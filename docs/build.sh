#!/bin/bash
# Build the Jekyll site
# Uses Jekyll 3.9.5 with sass-converter 1.5.2 to avoid ffi/sassc issues on arm64
# SCSS is pre-compiled to CSS via the `sass` CLI

set -e
cd "$(dirname "$0")"

# Recompile SCSS if source files changed
if [ -f _build.scss ]; then
  echo "Compiling SCSS..."
  sass --load-path _sass -I _sass --style compressed _build.scss assets/css/main.css
fi

echo "Building Jekyll site..."
# Temporarily hide Gemfile to prevent Bundler from interfering
mv Gemfile Gemfile.tmp 2>/dev/null || true
mv Gemfile.lock Gemfile.lock.tmp 2>/dev/null || true

# Default to 'build' if no subcommand provided
ARGS="${@:-build}"

# Output to 'docs' instead of default '_site'
DEST_ARGS="--destination docs"

ruby -e "
  gem 'jekyll', '4.4.1'
  gem 'jekyll-sass-converter'
  gem 'kramdown-parser-gfm', '1.1.0'
  gem 'rouge', '3.30.0'
  load Gem.bin_path('jekyll', 'jekyll', '4.4.1')
" -- $ARGS $DEST_ARGS
EXIT_CODE=$?

# Make site self-contained with relative URLs (only for build, not serve)
if [ "$EXIT_CODE" -eq 0 ] && [ "$ARGS" = "build" ] && [ -d "docs" ]; then
  ruby _scripts/make-relative.rb
fi

# Restore Gemfile
mv Gemfile.tmp Gemfile 2>/dev/null || true
mv Gemfile.lock.tmp Gemfile.lock 2>/dev/null || true

exit $EXIT_CODE