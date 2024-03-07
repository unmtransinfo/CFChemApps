// miscellaneous helper functions

// use compound name as download fname, if available
// fname will have been randomly generated (see utils.py)
function getDownloadName(name, fname) {
    var download_name = name;
    if (download_name === "---") {
        // --- == no name
        download_name = fname;
    }
    return download_name;
}