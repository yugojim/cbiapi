/*  
 * cbioEtl
*/
const express = require('express');
const router = express.Router();
const shell = require('shelljs');
const dt = require('date-and-time');
const fs = require('fs');

const apiName = 'cbioEtl';


// ETL and cbioPortal import
router.get('/import', (req, res) => {
  console.log('[' + apiName + '] Called, path=/');
  let vCaseID = "case-" + dt.format(new Date(), "YYYYMMDDHHmmss");

  shell.exec("sh cbioApiScript.sh " + vCaseID, (pStatus, pStdOut, pStdErr) => {
    if (pStdErr) console.error(pStdErr);
    console.log('[' + apiName + '] End, path=/');
  });

  res.json({
    caseId: vCaseID,
    status: "processing",
    message: "Please check after 30min."
  });
});

// Check the status of case ID
router.post('/status', (req, res) => {
  let vCaseID = req.body.caseId;
  if (fs.existsSync("./cbioApiStatus/" + vCaseID + ".1.running")) {
    if (fs.existsSync("./cbioApiStatus/" + vCaseID + ".2.etl")) {
      if (fs.existsSync("./cbioApiStatus/" + vCaseID + ".9.finish")) {
        res.json({ caseId: req.body.caseId, code: 9, status: "Case import completed." });
      } else {
        res.json({ caseId: req.body.caseId, code: 2, status: "cBioPortal data Importing." });
      }
    } else {
      res.json({ caseId: req.body.caseId, code: 1, status: "In ETL processing." });
    }
  } else {
    res.json({ caseId: req.body.caseId, code: -1, status: "Case not exists." });
  }
});

module.exports = router;