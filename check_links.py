#!/usr/bin/env python3
"""
Scan HTML files for broken links.
Checks if referenced files exist in the rangevoting.org directory.
"""

import os
import re
from pathlib import Path
from urllib.parse import urlparse, unquote

def extract_links(html_content):
    """Extract all href and src attributes from HTML content."""
    links = set()
    
    # Find all href attributes
    href_pattern = r'href=["\'](.*?)["\']'
    links.update(re.findall(href_pattern, html_content))
    
    # Find all src attributes
    src_pattern = r'src=["\'](.*?)["\']'
    links.update(re.findall(src_pattern, html_content))
    
    return links

def is_external_link(link):
    """Check if link is external (http, https, ftp, mailto, etc)."""
    return (
        link.startswith('http://') or
        link.startswith('https://') or
        link.startswith('ftp://') or
        link.startswith('mailto:') or
        link.startswith('//') or
        link.startswith('#') or
        link == ''
    )

def check_file_exists(link, base_dir, current_file_dir):
    """Check if a local file link exists."""
    # Parse the link (remove query strings and fragments)
    parsed = urlparse(link)
    path = unquote(parsed.path)
    
    if not path:
        return True  # Fragment-only links are always valid
    
    # Try absolute path first
    absolute_path = base_dir / path.lstrip('/')
    if absolute_path.exists():
        return True
    
    # Try relative to current file
    relative_path = current_file_dir / path
    if relative_path.exists():
        return True
    
    return False

def main():
    base_dir = Path('rangevoting.org')
    html_files = sorted(base_dir.glob('**/*.html'))
    
    print(f"Scanning {len(html_files)} HTML files for broken links...\n")
    
    broken_links = {}
    total_links = 0
    
    for html_file in html_files:
        try:
            with open(html_file, 'r', encoding='utf-8', errors='replace') as f:
                content = f.read()
            
            links = extract_links(content)
            total_links += len(links)
            
            file_broken = []
            
            for link in links:
                # Skip external links
                if is_external_link(link):
                    continue
                
                # Check if file exists
                if not check_file_exists(link, base_dir, html_file.parent):
                    file_broken.append(link)
            
            if file_broken:
                broken_links[str(html_file)] = file_broken
        
        except Exception as e:
            print(f"Error processing {html_file}: {e}")
    
    # Print results
    if broken_links:
        print(f"❌ Found {len(broken_links)} files with broken links:\n")
        
        for file_path, broken in sorted(broken_links.items()):
            print(f"  {file_path}")
            for link in sorted(set(broken)):
                print(f"    → {link}")
            print()
    else:
        print(f"✅ No broken links found!")
    
    print(f"\nSummary:")
    print(f"  Files scanned: {len(html_files)}")
    print(f"  Total links checked: {total_links}")
    print(f"  Files with broken links: {len(broken_links)}")

if __name__ == '__main__':
    main()
