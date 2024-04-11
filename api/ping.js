/*  
 * PING
*/
const express = require('express');
const router = express.Router();
const apiName = 'PING';

// INFO
router.get('/', (req, res) => {
    console.log('[' + apiName + '] Called, path=/');
    res.json({ greeting: 'Hello world' });
});


module.exports = router;