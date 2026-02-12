./gradlew :app:build
./gradlew :app:runSAR
# aws s3 sync app/out/ s3://${S3_BUCKET_NAME}/${S3_PATH}/