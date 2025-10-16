# RangeVoting.org Modernization Project

This repository contains the rangevoting.org website, currently being modernized for better maintainability and GitHub Pages deployment.

## Current Status

- ✅ Downloaded full site (955 HTML files)
- ✅ Extracted CSS from index.html to external stylesheet
- ✅ Set up local development environment with browser-sync
- 🔄 In progress: Batch processing remaining HTML files

## Local Development

### Starting the local server

```bash
./serve.sh
```

This will start browser-sync on `http://localhost:3000` with live reload enabled.

## Project Structure

```
rangevoting/
├── rangevoting.org/          # Original downloaded site
│   ├── index.html            # Modified to use external CSS
│   ├── assets/
│   │   └── css/
│   │       ├── original-layout.css  # Extracted original styles
│   │       └── style.css            # Modern clean styles (optional)
│   └── *.html                # 955+ other HTML files to process
├── serve.sh                  # Browser-sync development server
└── README.md                 # This file
```

## Strategy for Batch HTML Processing

Since there are 955 HTML files, we need a systematic approach:

### Phase 1: Analysis ✅

- [x] Extract inline CSS from index.html
- [x] Create external stylesheets
- [x] Set up development environment

### Phase 2: Batch CSS Extraction (Next)

We'll create a Python script to:

1. Find all HTML files with `<style>` tags
2. Extract inline CSS
3. Identify common patterns
4. Replace inline styles with `<link>` tags to external CSS

### Phase 3: CSS Consolidation

- Merge common styles across files
- Create a single master stylesheet
- Add CSS variables for easy theming
- Create alternate modern theme

### Phase 4: Cleanup & Enhancement

- Fix relative links
- Update deprecated HTML tags
- Add responsive meta tags
- Optimize images

## Batch Processing Tools

### Using Python + BeautifulSoup

```python
# Example script structure
from bs4 import BeautifulSoup
from pathlib import Path

def process_html_file(filepath):
    # Extract <style> blocks
    # Replace with <link> tags
    # Save modified HTML
    pass

# Process all files
for html_file in Path('rangevoting.org').glob('**/*.html'):
    process_html_file(html_file)
```

### Using sed/awk for bulk find-replace

```bash
# Example: Replace old links with new ones
find rangevoting.org -name "*.html" -exec sed -i '' 's/old-pattern/new-pattern/g' {} \;
```

### Using ast-grep for structural changes (if needed)

```bash
# For more complex HTML transformations
ast-grep --lang html -p '<style>' rangevoting.org/
```

## Deployment

The goal is to serve this as a static site on GitHub Pages:

1. Clean HTML files (no server-side code)
2. Relative paths for all assets
3. Single CSS file for consistency
4. All files in a `docs/` folder for GH Pages

## Next Steps

1. Create Python script to process remaining 954 HTML files
2. Extract and consolidate all CSS
3. Test across multiple pages
4. Create modern theme variant
5. Deploy to GitHub Pages

## Notes

- All 955 HTML files have inline `<style>` blocks
- The site uses old-school table layouts and floats
- Browser-sync watches for changes and auto-reloads
- Original content preserved, only styling modernized
