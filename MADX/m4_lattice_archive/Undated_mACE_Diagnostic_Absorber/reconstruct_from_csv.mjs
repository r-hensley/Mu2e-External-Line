#!/usr/bin/env node

/**
 * Reconstruct a runnable first-order MAD-X model from the received optics CSV.
 *
 * The CSV supplies beta, alpha, dispersion, and phase at the beginning,
 * midpoint, and end of each element.  Those quantities uniquely determine an
 * uncoupled 4x4 first-order transfer matrix for each row-to-row interval.
 * Dispersion
 * residuals determine the map's R16/R26/R36/R46 terms.  The generated MAD-X
 * model therefore reproduces the tabulated first-order optics without
 * pretending that the missing historical magnet definitions were recovered.
 */

import fs from "node:fs";
import path from "node:path";

const here = path.dirname(new URL(import.meta.url).pathname);
const csvPath = path.join(here, "mu2e-dabs-v2-functions.csv");
const mapDefinitionsPath = path.join(here, "mace_transport_maps.def");
const sequencePath = path.join(here, "mace_transport_reconstructed.seq");
const driverPath = path.join(here, "job_mace_transport_reconstructed.madx");
const parameterReportPath = path.join(here, "inferred_first_order_parameters.txt");

const numericColumns = new Set([
  "S[m]", "BX[m]", "BY[m]", "DX[m]", "DY[m]", "AX[1]", "AY[1]",
  "DDX[1]", "DDY[1]", "MUX[2*pi*rad]", "MUY[2*pi*rad]",
]);

/** Parse the simple comma-separated source table without altering it. */
function readCsv(filePath) {
  const lines = fs.readFileSync(filePath, "utf8").trim().split(/\r?\n/);
  const headers = lines.shift().split(",");
  return lines.map((line) => {
    const fields = line.split(",");
    return Object.fromEntries(headers.map((header, index) => [
      header,
      numericColumns.has(header) ? Number(fields[index]) : fields[index],
    ]));
  });
}

/** Return the 2x2 Courant-Snyder transport matrix between two table rows. */
function transverseMatrix(start, end, plane) {
  const betaKey = plane === "x" ? "BX[m]" : "BY[m]";
  const alphaKey = plane === "x" ? "AX[1]" : "AY[1]";
  const phaseKey = plane === "x" ? "MUX[2*pi*rad]" : "MUY[2*pi*rad]";
  const beta1 = start[betaKey];
  const beta2 = end[betaKey];
  const alpha1 = start[alphaKey];
  const alpha2 = end[alphaKey];
  const phase = 2 * Math.PI * (end[phaseKey] - start[phaseKey]);
  const cosine = Math.cos(phase);
  const sine = Math.sin(phase);

  return [
    Math.sqrt(beta2 / beta1) * (cosine + alpha1 * sine),
    Math.sqrt(beta1 * beta2) * sine,
    ((alpha1 - alpha2) * cosine - (1 + alpha1 * alpha2) * sine)
      / Math.sqrt(beta1 * beta2),
    Math.sqrt(beta1 / beta2) * (cosine - alpha2 * sine),
  ];
}

/** Build the uncoupled 4x4 matrix and its dispersion-driving column. */
function segmentMap(start, end) {
  const x = transverseMatrix(start, end, "x");
  const y = transverseMatrix(start, end, "y");
  const matrix = [
    [x[0], x[1], 0, 0],
    [x[2], x[3], 0, 0],
    [0, 0, y[0], y[1]],
    [0, 0, y[2], y[3]],
  ];

  const inputDispersion = [
    start["DX[m]"], start["DDX[1]"], start["DY[m]"], start["DDY[1]"],
  ];
  const outputDispersion = [
    end["DX[m]"], end["DDX[1]"], end["DY[m]"], end["DDY[1]"],
  ];
  const dispersionColumn = outputDispersion.map((value, row) => (
    value - matrix[row].reduce(
      (sum, coefficient, column) => sum + coefficient * inputDispersion[column],
      0,
    )
  ));

  // Complete the transverse-plus-dispersion map symplectically.  With
  // R = [[A,0,b],[c^T,1,r56],[0,0,1]], c = -A^T J b.
  const jTimesB = [
    dispersionColumn[1],
    -dispersionColumn[0],
    dispersionColumn[3],
    -dispersionColumn[2],
  ];
  const longitudinalRow = matrix[0].map((_, column) => -matrix.reduce(
    (sum, row, rowIndex) => sum + row[column] * jTimesB[rowIndex],
    0,
  ));

  return {
    length: end["S[m]"] - start["S[m]"],
    matrix,
    dispersionColumn,
    longitudinalRow,
  };
}

/** MAD-X-friendly scientific notation with enough digits for this inference. */
function fmt(value) {
  const clean = Math.abs(value) < 5e-17 ? 0 : value;
  return clean.toExponential(16).replace("e", "E");
}

/** Validate the strict begin/mid/end structure used by the reconstruction. */
function groupRows(rows) {
  if (rows.length % 3 !== 0) {
    throw new Error(`Expected row count divisible by three; found ${rows.length}`);
  }
  const groups = [];
  for (let index = 0; index < rows.length; index += 3) {
    const group = rows.slice(index, index + 3);
    const bases = group.map((row) => row.MARKER.replace(/_(BEG|MID|END)$/, ""));
    const suffixes = group.map((row) => row.MARKER.match(/_(BEG|MID|END)$/)?.[1]);
    if (new Set(bases).size !== 1 || suffixes.join(",") !== "BEG,MID,END") {
      throw new Error(`Malformed row triplet beginning at CSV row ${index + 2}`);
    }
    if (new Set(group.map((row) => row.CATEGORY)).size !== 1) {
      throw new Error(`Category changes inside ${bases[0]}`);
    }
    groups.push(group);
  }
  for (let index = 1; index < groups.length; index += 1) {
    const previousEnd = groups[index - 1][2]["S[m]"];
    const currentBegin = groups[index][0]["S[m]"];
    if (Math.abs(previousEnd - currentBegin) > 5e-10) {
      throw new Error(`Longitudinal discontinuity before ${groups[index][0].MARKER}`);
    }
  }
  return groups;
}

/** Derive a hard-edge K1 estimate from a tabulated quadrupole's full map. */
function inferredQuadrupoleK1(group) {
  const full = segmentMap(group[0], group[2]);
  const x = full.matrix;
  const kx = -x[1][0] / x[0][1];
  const ky = -x[3][2] / x[2][3];
  return { k1: (kx - ky) / 2, xEstimate: kx, yEstimate: -ky };
}

/** Estimate a sector-bend kick vector from the dispersion residual. */
function inferredDipole(group) {
  const full = segmentMap(group[0], group[2]);
  const kickX = full.dispersionColumn[1];
  const kickY = full.dispersionColumn[3];
  const sineMagnitude = Math.min(1, Math.hypot(kickX, kickY));
  const dominantKick = Math.abs(kickX) >= Math.abs(kickY) ? kickX : kickY;
  const signedAngle = Math.sign(dominantKick || 1) * Math.asin(sineMagnitude);
  const sign = Math.sign(signedAngle || 1);
  const tilt = Math.atan2(kickY / sign, kickX / sign);
  return { angle: signedAngle, tilt, kickX, kickY };
}

const rows = readCsv(csvPath);
const groups = groupRows(rows);
const mapDefinitions = [
  "! Generated 2026-08-15 by reconstruct_from_csv.mjs.",
  "! First-order maps inferred from the rounded CSV optics samples.",
  "CSV_POINT: marker;",
  "",
];
const sequenceLines = [
  "! Generated 2026-08-15 by reconstruct_from_csv.mjs.",
  "! Marker names and coordinates follow the received CSV exactly.",
  `MACE_DABS: sequence, l=${fmt(rows.at(-1)["S[m]"])}, refer=centre;`,
];
const reportLines = [
  "INFERRED FIRST-ORDER PARAMETERS",
  "Generated 2026-08-15 from mu2e-dabs-v2-functions.csv",
  "",
  "These are estimates derived from rounded Twiss/phase/dispersion samples.",
  "They are not recovered historical settings or control-system values.",
  "",
  "QUADRUPOLES",
  "name                         length_m        inferred_K1_m^-2   x_plane_estimate   y_plane_estimate",
];

let mapNumber = 0;
sequenceLines.push(`  ${rows[0].MARKER}: CSV_POINT, at=${fmt(rows[0]["S[m]"])};`);
for (let index = 0; index < rows.length - 1; index += 1) {
  const start = rows[index];
  const finish = rows[index + 1];
  const map = segmentMap(start, finish);
  mapNumber += 1;
  const className = `MAP_${String(mapNumber).padStart(4, "0")}`;
  const instanceName = `E_${String(mapNumber).padStart(4, "0")}`;
  const m = map.matrix;
  const b = map.dispersionColumn;
  const c = map.longitudinalRow;
  mapDefinitions.push(
    `${className}: matrix, l=${fmt(map.length)},`,
    `  rm11=${fmt(m[0][0])}, rm12=${fmt(m[0][1])}, rm16=${fmt(b[0])},`,
    `  rm21=${fmt(m[1][0])}, rm22=${fmt(m[1][1])}, rm26=${fmt(b[1])},`,
    `  rm33=${fmt(m[2][2])}, rm34=${fmt(m[2][3])}, rm36=${fmt(b[2])},`,
    `  rm43=${fmt(m[3][2])}, rm44=${fmt(m[3][3])}, rm46=${fmt(b[3])},`,
    `  rm51=${fmt(c[0])}, rm52=${fmt(c[1])}, rm53=${fmt(c[2])}, rm54=${fmt(c[3])},`,
    "  rm55=1, rm66=1;",
    "",
  );
  const center = (start["S[m]"] + finish["S[m]"]) / 2;
  sequenceLines.push(`  ${instanceName}: ${className}, at=${fmt(center)};`);
  sequenceLines.push(`  ${finish.MARKER}: CSV_POINT, at=${fmt(finish["S[m]"])};`);
}

for (const group of groups) {
  const [begin, , end] = group;
  if (begin.CATEGORY === "QUADRUPOLE") {
    const estimate = inferredQuadrupoleK1(group);
    reportLines.push(
      `${begin.MARKER.replace(/_BEG$/, "").padEnd(28)} `
      + `${(end["S[m]"] - begin["S[m]"]).toFixed(6).padStart(12)} `
      + `${estimate.k1.toFixed(10).padStart(21)} `
      + `${estimate.xEstimate.toFixed(10).padStart(18)} `
      + `${estimate.yEstimate.toFixed(10).padStart(18)}`,
    );
  }
}
sequenceLines.push("endsequence;", "");

reportLines.push(
  "",
  "DIPOLES",
  "name                         length_m        inferred_angle_rad  inferred_tilt_rad   R26_x             R46_y",
);
for (const group of groups.filter((group) => group[0].CATEGORY === "DIPOLE")) {
  const [begin, , end] = group;
  const estimate = inferredDipole(group);
  reportLines.push(
    `${begin.MARKER.replace(/_BEG$/, "").padEnd(28)} `
    + `${(end["S[m]"] - begin["S[m]"]).toFixed(6).padStart(12)} `
    + `${estimate.angle.toFixed(10).padStart(21)} `
    + `${estimate.tilt.toFixed(10).padStart(18)} `
    + `${estimate.kickX.toExponential(7).padStart(17)} `
    + `${estimate.kickY.toExponential(7).padStart(17)}`,
  );
}

const first = rows[0];
const driver = `! Local first-order reconstruction generated 2026-08-15.
! This is not recovered historical source.  See !NOTE.md.
beam, particle=proton, energy=8.93828;

call, file="mace_transport_maps.def";
call, file="mace_transport_reconstructed.seq";

csv_start: beta0,
  betx=${fmt(first["BX[m]"])}, alfx=${fmt(first["AX[1]"])},
  bety=${fmt(first["BY[m]"])}, alfy=${fmt(first["AY[1]"])},
  dx=${fmt(first["DX[m]"])}, dpx=${fmt(first["DDX[1]"])},
  dy=${fmt(first["DY[m]"])}, dpy=${fmt(first["DDY[1]"])};

use, sequence=MACE_DABS;
set, format="20.12g";
select, flag=twiss, clear;
select, flag=twiss,
  column=name,s,betx,alfx,bety,alfy,dx,dpx,dy,dpy,mux,muy,keyword,l;
twiss, beta0=csv_start, sequence=MACE_DABS,
  file="twiss_mace_transport_reconstructed.tfs";
stop;
`;

fs.writeFileSync(mapDefinitionsPath, `${mapDefinitions.join("\n")}\n`);
fs.writeFileSync(sequencePath, `${sequenceLines.join("\n")}\n`);
fs.writeFileSync(driverPath, driver);
fs.writeFileSync(parameterReportPath, `${reportLines.join("\n")}\n`);

console.log(`Validated ${rows.length} CSV rows in ${groups.length} triplets.`);
console.log(`Generated ${mapNumber} row-to-row transfer maps.`);
console.log(`Wrote ${path.basename(mapDefinitionsPath)}.`);
console.log(`Wrote ${path.basename(sequencePath)}.`);
console.log(`Wrote ${path.basename(driverPath)}.`);
console.log(`Wrote ${path.basename(parameterReportPath)}.`);
