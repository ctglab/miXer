#!/bin/bash

# Check if the keyword is provided
if [ -z "$1" ]; then
    echo "Please provide a keyword to select operating mode."
    exit 1
fi

# Check if the JSON file is provided
if [ -z "$2" ]; then
    echo "Please provide a JSON file with run parameters."
    exit 1
fi

KEYWORD="$1"
JSON_FILE="$2"


# Check if the JSON file exists
if [ ! -f "$JSON_FILE" ]; then
    echo "JSON file not found: $JSON_FILE"
    exit 1
fi

# Check the keyword and run the corresponding script with the JSON file
case $KEYWORD in
    "script1")
        exec /entrypoints/mixer_full_inference.sh "$JSON_FILE"
        ;;
    "script2")
        exec /entrypoints/mixer_preprocessing_and_inference.sh "$JSON_FILE"
        ;;
    "script3")
        exec /entrypoints/mixer_preprocessing.sh "$JSON_FILE"
        ;;
    *)
        echo "Invalid keyword: $KEYWORD"
        exit 1
        ;;
esac
