#!/bin/bash

# This assumes that you're in the "build" of a cmake build process. E.g.,
# you've already done something like:
#
#   mkdir build
#   cd build
#   cmake ..
#
# Then this will (try to) circumvent cmake's aggressive (i.e., predefined macro
# oblivious) dependency scanner by calling makedepend directly for every .cpp
# file of every target.

NUM_THREADS="1"
while getopts ":C:j:h" opt; do
  case $opt in
    C)
      if ! cd "$OPTARG" 2>/dev/null
      then
        (>&2 echo "Failed to change directory to $OPTARG")
        exit 1
      fi
      ;;
    h)
      echo "
Usage:
  
    make.sh [-j #] [-C dir] [project]"
      exit 1
      ;;
    j)
      NUM_THREADS="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done

# Shift so that $# makes sense
shift $((OPTIND-1))

if [ ! -f "Makefile" ]; then
  (>&2 echo "Makefile not found")
  exit 1
fi

# I don't want to do this if I don't have to:
#make depend

if [[ $# -eq 0 ]] ; then
  # Try to get project name
  PROJECTS=`find . -type d -depth 1 | sed -e "s/.\///g" | grep -v CMakeFiles | grep -v "^\."`
else
  PROJECTS=`echo "$@" | tr ' ' '\n'`
fi

for PROJECT in $PROJECTS
do

  PROJECT="${PROJECT%/}"
  
  if [ ! -d "$PROJECT" ]; then
    (>&2 echo "$PROJECT directory not found")
    exit 1
  fi
  
  
  for TARGET_DIR in `find "$PROJECT/CMakeFiles" -depth 1 | grep "\.dir"`
  do
    TARGET_DIR="${TARGET_DIR%/}"
    if ! grep --quiet "# make\.sh was here" "$TARGET_DIR/depend.make"
    then
      if grep --quiet "Empty dependencies" "$TARGET_DIR/depend.make"
      then
        make -f "$TARGET_DIR/build.make" "$TARGET_DIR/depend"
      fi
      if grep --quiet "Empty dependencies" "$TARGET_DIR/depend.make"
      then
        (>&2 echo "cmake's make failed to build a dependency file")
        (>&2 echo "[Consider issuing \`make depend\`]")
        continue
      fi
    fi
    TARGET=`basename -s .dir $TARGET_DIR`
    echo -e "\033[1;35mHacking dependencies of target $TARGET\033[0m"
    SRC_FILES=`cat "$TARGET_DIR/depend.make" | sed -n "s~^#*\($TARGET_DIR\)\(.*\)\.o: \2~\2~gp"`
    if grep --quiet "CMAKE generated file" "$TARGET_DIR/depend.make"
    then
      mv "$TARGET_DIR/depend.make"{,.bk}
    fi
    # Create new, empty depend.make file
    echo "# make.sh was here" > "$TARGET_DIR/depend.make"
    for SRC in $SRC_FILES
    do
      EXT="${SRC##*.}"
      makedepend -DIGL_STATIC_LIBRARY -p "$TARGET_DIR" -o .${EXT}.o -a -f "$TARGET_DIR/depend.make" -w0 $SRC 2>/dev/null 
      # Add commented self so that next run of make.sh can use this depend.make file
      echo "#$TARGET_DIR$SRC.o: $SRC" >> "$TARGET_DIR/depend.make"
    done
    #echo -e "\033[0;35mHanding off to make $TARGET\033[0m"
    make -j$NUM_THREADS -f "${TARGET_DIR}/build.make" "${TARGET_DIR}/build" 
  done
done
