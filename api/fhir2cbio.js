/*  
 * FHIR ETL to cBio
*/
const express = require('express');
const router = express.Router();
const shell = require('shelljs');
const dt = require('date-and-time');
const fs = require('fs');
const CONF = require('../cbioApiConfig');

const apiName = 'fhir2cbio';

// Trigger when user uploaded report in FHIR portal
router.post('/up', (req, res) => {
  console.log('[' + apiName + '] Called, path=/up');
  let vCaseID = "case-" + dt.format(new Date(), "YYYYMMDDHHmmss");

  // // start saving csv
  
  // Formulate the path for saving the CSV file
  const saveCsvFolder = 'save_csv';
  const csvPath = `${saveCsvFolder}/${vCaseID}.ngs.csv`
  console.log('Saving NGS.csv file to path:', csvPath);

  // Save the CSV file
  // fs.writeFileSync(csvPath, req.body.NGS_csv, 'utf-8');
  fs.writeFileSync(csvPath, req.body.NGS_csv);
  console.log('NGS.csv file saved successfully.');

  // // finish saving csv

  shell.exec("sh ./sh/testing_" + CONF.environment + ".sh " + vCaseID, (pStatus, pStdOut, pStdErr) => {
    if (pStdErr) console.error(pStdErr);
    console.log('[' + apiName + '] End, path=/up');
  });

  res.json({
    caseId: vCaseID,
    status: "processing",
    message: "Please check after 10sec.",
    arg: req.body
  });
});



// Check the status of case ID
router.post('/status', (req, res) => {
  let vCaseID = req.body.caseId;
  if (fs.existsSync( CONF.path_docker + "/cbioApiStatus/" + vCaseID + ".1.ready_to_start")) {
    if (fs.existsSync(CONF.path_docker + "/cbioApiStatus/" + vCaseID + ".2.finish_running_etl_process")) {
      if (fs.existsSync(CONF.path_docker + "/cbioApiStatus/" + vCaseID + ".9.finish_uploading")) {
        res.json({ caseId: req.body.caseId, code: 9, status: "New cases import completed." });
      } else {
        res.json({ caseId: req.body.caseId, code: 2, status: "Data to cBioPortal is still importing." });
      }
    } else {
      res.json({ caseId: req.body.caseId, code: 1, status: "Data is still in ETL processing." });
    }
  } else {
    res.json({ caseId: req.body.caseId, code: -1, status: "Case not exists." });
  }
});


// Trigger when user uploaded report in FHIR portal
router.post('/addToCbio', (req, res) => {
  console.log('[' + apiName + '] Called, path=/addToCbio');
  let vCaseID = "case-" + dt.format(new Date(), "YYYYMMDDHHmmss");

  // // start saving csv
  
  // Formulate the path for saving the text file
  const saveCsvFolder = 'save_csv';
  const csvPath = `${saveCsvFolder}/${vCaseID}.ngs.txt`
  console.log('Saving NGS.txt file to path:', csvPath);

  // Save the CSV file
  // fs.writeFileSync(csvPath, req.body.NGS_csv, 'utf-8');
  fs.writeFileSync(csvPath, JSON.stringify(req.body));
  console.log('NGS.txt file saved successfully.');

  // finish saving txt
  
  // shell.exec("sh ./sh/testing_" + CONF.environment + ".sh " + vCaseID, (pStatus, pStdOut, pStdErr) => {
  //   if (pStdErr) console.error(pStdErr);
  //   console.log('[' + apiName + '] End, path=/up');
  // });

  res.json({
    caseId: vCaseID,
    status: "processing",
    message: "Please check after 10sec.",
    arg: req.body
  });
});


module.exports = router;