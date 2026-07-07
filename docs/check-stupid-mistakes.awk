#! /usr/bin/awk -f

# /^```@/ {
#     print
# }

BEGIN {
    valid_block_tags["contents"] = 0
    valid_block_tags["meta"] = 0
    valid_block_tags["repl"] = 0
    valid_block_tags["docs"] = 0
}

match($0, /^```@(.+)$/, ms) {
    m = ms[1]
    if (!(m in valid_block_tags))
        print "invalid tag: " m " --- (" FILENAME ":" NR ")"
}
