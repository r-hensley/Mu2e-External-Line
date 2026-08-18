#!/usr/bin/env node

/** Compare the received CSV sample markers with a regenerated MAD-X TFS. */

import fs from "node:fs";
import path from "node:path";

const here = path.dirname(new URL(import.meta.url).pathname);
const csvPath = path.join(here, "mu2e-dabs-v2-functions.csv");
const tfsPath = path.resolve(process.argv[2] ?? path.join(
  here,
  "twiss_mace_transport_reconstructed.tfs",
));

/** Read the source CSV and convert numerical fields to numbers. */
function readCsv(filePath) {
  const lines = fs.readFileSync(filePath, "utf8").trim().split(/\r?\n/);
  const headers = lines.shift().split(",");
  return lines.map((line) => {
    const fields = line.split(",");
    return Object.fromEntries(headers.map((header, index) => [
      header,
      index < 2 ? fields[index] : Number(fields[index]),
    ]));
  });
}

/** Read the whitespace-delimited MAD-X TFS data table. */
function readTfs(filePath) {
  const lines = fs.readFileSync(filePath, "utf8").split(/\r?\n/);
  let columns = [];
  const rows = [];
  for (const line of lines) {
    if (line.startsWith("*")) {
      columns = line.slice(1).trim().split(/\s+/);
      continue;
    }
    if (!columns.length || !line.trim() || line.startsWith("@") || line.startsWith("$")) {
      continue;
    }
    const fields = line.trim().match(/"[^"]*"|\S+/g);
    if (!fields || fields.length < columns.length) {
      continue;
    }
    rows.push(Object.fromEntries(columns.map((column, index) => [
      column,
      fields[index].replace(/^"|"$/g, ""),
    ])));
  }
  return rows;
}

const columnMap = new Map([
  ["S[m]", "S"],
  ["BX[m]", "BETX"],
  ["AX[1]", "ALFX"],
  ["BY[m]", "BETY"],
  ["AY[1]", "ALFY"],
  ["DX[m]", "DX"],
  ["DDX[1]", "DPX"],
  ["DY[m]", "DY"],
  ["DDY[1]", "DPY"],
  ["MUX[2*pi*rad]", "MUX"],
  ["MUY[2*pi*rad]", "MUY"],
]);

const csvRows = readCsv(csvPath);
const tfsRows = readTfs(tfsPath);
const tfsByName = new Map(tfsRows.map((row) => [row.NAME.toUpperCase(), row]));
const missing = csvRows.filter((row) => !tfsByName.has(row.MARKER.toUpperCase()));

console.log("MACE CSV TO MAD-X RECONSTRUCTION COMPARISON");
console.log(`Generated: 2026-08-15`);
console.log(`CSV: ${path.basename(csvPath)}`);
console.log(`TFS: ${path.basename(tfsPath)}`);
console.log("");
console.log(`CSV sample rows: ${csvRows.length}`);
console.log(`MAD-X TFS rows: ${tfsRows.length}`);
console.log(`Matched named CSV samples: ${csvRows.length - missing.length}`);
console.log(`Missing named CSV samples: ${missing.length}`);
console.log("");
console.log("Maximum absolute differences at matched sample markers:");
console.log("CSV column              TFS column    maximum difference       marker");

for (const [csvColumn, tfsColumn] of columnMap) {
  let maximum = -1;
  let maximumMarker = "";
  for (const csvRow of csvRows) {
    const tfsRow = tfsByName.get(csvRow.MARKER.toUpperCase());
    if (!tfsRow) continue;
    const difference = Math.abs(csvRow[csvColumn] - Number(tfsRow[tfsColumn]));
    if (difference > maximum) {
      maximum = difference;
      maximumMarker = csvRow.MARKER;
    }
  }
  console.log(
    `${csvColumn.padEnd(23)} ${tfsColumn.padEnd(12)} `
    + `${maximum.toExponential(12).padStart(20)}    ${maximumMarker}`,
  );
}

const allMatchAtSourcePrecision = [...columnMap].every(([csvColumn, tfsColumn]) => (
  csvRows.every((csvRow) => {
    const tfsRow = tfsByName.get(csvRow.MARKER.toUpperCase());
    return tfsRow && Number(tfsRow[tfsColumn]).toFixed(6) === csvRow[csvColumn].toFixed(6);
  })
));

console.log("");
console.log(`All values match at the CSV's six-decimal precision: ${allMatchAtSourcePrecision}`);
if (missing.length) {
  console.log(`Missing markers: ${missing.map((row) => row.MARKER).join(", ")}`);
  process.exitCode = 1;
}
if (!allMatchAtSourcePrecision) process.exitCode = 1;
