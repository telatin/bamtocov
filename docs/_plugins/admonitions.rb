require 'kramdown'

Jekyll::Hooks.register :documents, :pre_render do |doc|
  # Convert Docusaurus-style admonitions to HTML
  # Supports: :::note, :::tip, :::caution, :::info, :::warning
  # Also supports titled variants like :::info[Title]
  doc.content = doc.content.gsub(/^:::(note|tip|caution|info|warning)(?:\[([^\]]*)\])?\s*\n(.*?)\n:::/m) do
    type = $1
    title = $2 || type.capitalize
    content = $3.strip

    # Process markdown content through kramdown
    content_html = Kramdown::Document.new(content).to_html.strip

    icon = case type
           when 'note' then '📝'
           when 'tip' then '💡'
           when 'caution', 'warning' then '⚠️'
           when 'info' then 'ℹ️'
           else '📝'
           end

    <<~HTML
      <div class="admonition admonition-#{type}">
        <div class="admonition-heading">
          <span class="admonition-icon">#{icon}</span>
          <span class="admonition-title">#{title}</span>
        </div>
        <div class="admonition-content">
      #{content_html}
        </div>
      </div>
    HTML
  end
end
