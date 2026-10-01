#!/bin/bash

# TODO: Provide an option to compare exported results with Mathematica and FORM.

KIRA="$1"
KIRA_JOB_ID="$2"
SYMCOMPARE="$3"
INPUT_PROJCT_DIR="$4/$6"
REFERENCE_PROJECT_DIR="$5/$6"
TEST_DIR="`pwd`/tests"
TEST_PROJECT_DIR="$TEST_DIR/$6"
TEST_PROJECT_TAG=".kira_test_project_tag"
TEST_LEVEL="$7"

if [ ! -x "$KIRA" ]; then
  echo "ERROR: Kira path does not point to an executable file."
  exit 1
fi
if [ -z "$KIRA_JOB_ID" ]; then
  echo "ERROR: Kira job file missing."
  exit 1
fi
if [ ! -x "$SYMCOMPARE" ]; then
  echo "ERROR: symcompare path does not point to an executable file."
  exit 1
fi
if [ ! -d "$INPUT_PROJCT_DIR" ]; then
  echo "ERROR: Input project directory doesn't exist."
  exit 1
fi
if [ ! -d "$REFERENCE_PROJECT_DIR" ]; then
  echo "ERROR: Reference project directory doesn't exist."
  exit 1
fi
if [ ! -d "$TEST_DIR" ]; then
  echo "ERROR: Test directory doesn't exist."
  exit 1
fi
if [ -z "$TEST_LEVEL" ]; then
  TEST_LEVEL="3"
fi

if [ -f "$TEST_PROJECT_DIR/$TEST_PROJECT_TAG" ]; then
  echo "Delete old test project directory $TEST_PROJECT_DIR"
  rm -r "$TEST_PROJECT_DIR"
elif [ -e "$TEST_PROJECT_DIR" ]; then
  echo "$TEST_PROJECT_DIR exists, but is not a test project directory."
  echo "I better quit before I delete something valuable."
  exit 1
fi
cp -r "$INPUT_PROJCT_DIR" "$TEST_DIR/"
touch "$TEST_PROJECT_DIR/$TEST_PROJECT_TAG"

cd "$TEST_PROJECT_DIR"
"$KIRA" "job_${KIRA_JOB_ID}_${TEST_LEVEL}.yaml"
"$SYMCOMPARE" "--mode=$TEST_LEVEL" "$TEST_PROJECT_DIR" "$REFERENCE_PROJECT_DIR" "$KIRA_JOB_ID"
