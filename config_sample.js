
// cp this config to ./cbioApiConfig.js
// environmentName range in { 'uat', 'prod', 'dev' }
let environmentName = 'dev';

module.exports = require('./apiConfig/conf_' + environmentName);
