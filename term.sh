#!/bin/bash
read -r -d '' script <<'EOF'
on run argv
tell application "iTerm"
    activate
    set myterm to (create window with default profile)
    tell myterm
        tell the current session
            -- write text argv
            -- repeat with arg in argv
               -- say arg
               write text "cd /Users/natj/projects/radpic"
               write text "gdb radpic"
            -- end repeat
        end tell
    end tell
end tell
end run
EOF
echo "$script" | osascript ``-'' $@
