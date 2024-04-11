const app = require('express')();
const bodyParser = require('body-parser');
const CONF = require('./cbioApiConfig');
const cors = require('cors')

const fs = require('fs');
const https = require('https');

app.use(cors({
  origin: '*'
}));

app.use(bodyParser.json());

app.use('/ping', require('./api/ping'));
// app.use('/i/cbio', require('./cbioApiETL'));
app.use('/fhir', require('./api/fhir2cbio'));

if (CONF.environment == 'PROD') {
  https.createServer({
    key: fs.readFileSync('./ssl/xxx.com_private_key.key'),
    cert: fs.readFileSync('./ssl/xxx.com_ssl_certificate.pem'),
  }, app).listen(CONF.servicePortSSL, startListenSSL);
} else {
  app.listen(CONF.servicePort, startListen);
}

function startListen() {
  console.log('API launched: Now listening on port(' + CONF.servicePort + ').');
}

function startListenSSL() {
  console.log('API launched: Now listening with SSL on port(' + CONF.servicePortSSL + ').');
}
