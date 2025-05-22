RED='\033[1;91m'
SKY='\033[1;96m'
PURPLE='\033[1;95m'
Color_reset='\033[0m'

if [ $# -ne 1 ]; then
  echo -e "${RED}LISTING PRESETS ${SKY}"
  # cmake --list-presets
  cmake --build --list-presets
  echo -e "${RED}TO BUILD, DO: ${PURPLE}bash $0 <preset>${Color_reset}"
  exit 1
fi

cmake --preset $1
cmake --build --preset $1
