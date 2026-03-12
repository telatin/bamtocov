#!/usr/bin/env ruby
# Convert absolute URLs to relative URLs for self-contained local browsing
# Run after Jekyll build: ruby _scripts/make-relative.rb

require 'fileutils'

SITE_DIR = 'docs'

def depth_prefix(file_path)
  # Calculate relative path prefix based on file depth from _site root
  relative = file_path.sub(%r{^#{SITE_DIR}/?}, '')
  depth = relative.count('/')
  depth == 0 ? './' : ('../' * depth)
end

def process_html(file_path)
  content = File.read(file_path, encoding: 'UTF-8')
  prefix = depth_prefix(file_path)

  # Replace absolute paths with relative ones
  # Match href="/...", src="/...", url('/...'), url("/...")
  modified = content.gsub(%r{(href|src|content)=["']/(?!/)}, "\1=\"#{prefix}")
  modified = modified.gsub(%r{url\(["']?/(?!/)}i, "url(#{prefix}")

  # Handle meta refresh: content="0; url=/docs/" -> content="0; url=./docs/"
  modified = modified.gsub(%r{(content=["'][^"']*url=)/(?!/)}i, "\1#{prefix}")

  # Fix any double slashes that might occur (except http://)
  modified = modified.gsub(%r{(?<!:)//+}, '/')

  if content != modified
    File.write(file_path, modified)
    puts "  Fixed: #{file_path}"
  end
end

def process_css(file_path)
  content = File.read(file_path, encoding: 'UTF-8')
  prefix = depth_prefix(file_path)

  # Replace url(/...) patterns in CSS
  modified = content.gsub(%r{url\(["']?/(?!/)}i, "url(#{prefix}")


  modified = modified.gsub(%r{(?<!:)//+}, '/')

  if content != modified
    File.write(file_path, modified)
    puts "  Fixed: #{file_path}"
  end
end

puts "Making site self-contained (relative URLs)..."

Dir.glob("#{SITE_DIR}/**/*.html").each { |f| process_html(f) }
Dir.glob("#{SITE_DIR}/**/*.css").each { |f| process_css(f) }

puts "Done!"
