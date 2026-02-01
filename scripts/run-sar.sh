#!/bin/bash

# .env ファイルから設定を読み込む
if [ -f .env ]; then
    export $(cat .env | grep -v '^#' | xargs)
else
    echo "Error: .env file not found. Please create .env file based on .env.example"
    exit 1
fi

# 必須変数のチェック
if [ -z "$S3_BUCKET_NAME" ] || [ -z "$S3_PATH" ]; then
    echo "Error: S3_BUCKET_NAME and S3_PATH must be set in .env file"
    exit 1
fi

./gradlew :app:build
./gradlew :app:runSAR
aws s3 sync app/out/ s3://${S3_BUCKET_NAME}/${S3_PATH}/