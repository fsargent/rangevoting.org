#!/bin/bash

# Simple browser-sync server for local development
# Watches for changes and auto-reloads the browser

echo "Starting browser-sync server for rangevoting.org..."
echo "Server will be available at http://localhost:3000"
echo "Press Ctrl+C to stop"
echo ""

cd rangevoting.org

browser-sync start \
	--server \
	--files "**/*.html, **/*.css, **/*.js, **/*.png, **/*.jpg, **/*.gif" \
	--port 3000 \
	--no-notify \
	--no-open
