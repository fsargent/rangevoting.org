# RangeVoting.org

A comprehensive resource for understanding range voting and voting system reform.

## About

**RangeVoting.org** is dedicated to promoting range voting as a superior alternative to traditional voting systems. Range voting (also called score voting) allows voters to score each candidate on a scale (e.g., 0-9), enabling more nuanced expression of voter preferences than single-choice systems.

Founded and maintained by voting system researchers, this site provides:
- Educational content on voting theory
- Comparative analysis of different voting systems
- Research papers and technical publications
- Interactive tools and calculators
- Arguments for voting reform

## Site Structure

### Core Pages
- **index.html** - Main page with navigation and introduction to range voting
- **RangeVoting.html** - Comprehensive range voting resource
- **IEVS.html** - Institute for Election Voting Systems information
- **VermontRept.html** - Vermont election-related reports

### Content Categories

#### Voting Systems
- Comparison of different voting methods (IRV, Borda, Approval, etc.)
- Analysis of voting system properties
- Mathematical foundations of voting theory

#### Research & Publications
- **WarrenSmithPages/** - Warren D. Smith's academic works
  - works.html - Comprehensive bibliography of research papers
  - Technical papers on voting systems, physics, mathematics, and computer science

#### Tools
- **VoteCostCalc.html** - Interactive voting machine cost calculator
  - Compare manual vs. machine vote counting costs
  - Contributed by Brian Olson

#### Asset Organization
- **assets/css/** - Stylesheets
  - original-layout.css - Main site styling
- **assets/documents/** - PDFs, PostScript files, and text documents
- **assets/images/** - PNG, JPG, GIF images
- **assets/archive/** - Archived HTML pages no longer actively linked

## Technical Details

### Architecture
- Static HTML site (no server-side processing required)
- Clean CSS styling with responsive design
- Client-side JavaScript for interactive tools
- Organized asset structure for easy maintenance

### Key Features
- **Navigation Sidebar** - Table of contents for easy browsing
- **Responsive Layout** - Works on desktop and mobile devices
- **Professional Typography** - Clean, readable design
- **Interactive Calculator** - Real-time voting cost analysis

### Recent Improvements
- Modernized HTML structure (HTML5)
- Improved CSS organization
- Refactored legacy CGI tools to client-side JavaScript
- Organized media assets into appropriate directories
- Enhanced navigation and content hierarchy

## Asset Organization

### Media Files
All media files are organized in the `assets/` directory:

```
assets/
├── css/                 # Stylesheets
├── documents/          # PDFs, PostScript, text files
├── images/             # PNG, JPG, GIF files
├── archive/            # Archived/unused content
├── audio/              # Audio files
├── data/               # Data files
├── graphics/           # Additional graphics
├── presentations/      # Presentation files
├── source/             # Source code
└── spreadsheets/       # Data spreadsheets
```

### Root Directory
The root directory contains only HTML files and subdirectories, keeping the site organized and maintainable.

## Contributing

### Adding New Content
1. Create new HTML files in the root directory
2. Link to them from the navigation (index.html)
3. Place any associated assets in the appropriate `assets/` subdirectory
4. Update links to use the new asset paths

### Updating Styles
- Edit `assets/css/original-layout.css` for global styling
- Keep inline styles to a minimum
- Use semantic HTML5 elements

### Adding Pages
- Use HTML5 DOCTYPE and semantic tags
- Include the site's CSS stylesheet
- Follow the existing navigation pattern
- Ensure all links are relative and point to correct asset locations

## Link Conventions

### Internal Links
- Root files: `<a href="PageName.html">`
- Assets: `<a href="/assets/documents/file.pdf">`
- Warren Smith papers: `<a href="WarrenSmithPages/homepage/works.html">`

### Relative Paths
All links should be relative to avoid breaking when site is served from different domains.

## Tools & Resources

### Interactive Tools
- **Voting Cost Calculator** (VoteCostCalc.html) - Analyzes manual vs. machine counting costs

### External Resources
- Brian Olson's voting reform site: http://bolson.org/voting/
- Warren D. Smith's academic page: http://math.temple.edu/~wds/homepage/

## Maintenance

### Regular Tasks
- Check for broken links periodically
- Update outdated information in papers/resources
- Maintain asset organization as new content is added

### File Naming
- Use descriptive, CamelCase names for HTML files
- Avoid query parameters in filenames
- Keep filenames concise but informative

## License & Attribution

Content on this site includes:
- Original research by voting reform advocates
- Academic papers by Warren D. Smith and contributors
- Tools created by Brian Olson and others
- Public domain and licensed content

Please see individual papers and resources for specific licensing information.

## Contact & Support

For information about range voting and voting system reform, visit the main site or contact the site administrators.

---

**Last Updated:** October 2025

This README documents the structure and organization of RangeVoting.org to help maintain and expand the site.
