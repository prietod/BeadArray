# join <delim> <list...>
# Concatenates the list elements with the delimiter passed as first parameter
#
# Ex: join , a b c
#  -> a,b,c
function join {
  local IFS="$1"
  shift
  echo "$*"
}

function render-template {
  eval "echo \"$(cat $1)\""
}


