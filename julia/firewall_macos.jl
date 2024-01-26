# On a Mac with firewall enabled, running MAGEMin will result 
# in a popup window that says:
#       "Accept incoming connections" which yo should Allow or Deny
#
# This is a bit annoying, so this julia script fixes that
#
# Note that you ust have administrator rights on your machine as we need to run "sudo"
#
# Run this script from the terminal with
#
# julia firewall_macos.jl  

using MAGEMin_jll

firewall_app = "/usr/libexec/ApplicationFirewall/socketfilterfw"

# 1) Deactivate firewall
run(`sudo $firewall_app --setglobalstate off`) 

# 2) Add MAGEMin executable to firewall 
exe = MAGEMin_jll.MAGEMin_path
run(`sudo $firewall_app --add $(exe)`) 

# 3) Block incoming connections
run(`sudo $firewall_app --block $(exe)`) 

# 4) Activate firewall again
run(`sudo $firewall_app --setglobalstate on`) 

