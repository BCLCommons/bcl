
def GetHTMLHeader():
  return """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">

<head>
<script src="sorttable.js"></script>
</head>

<html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
"""

def GetHTMLFooter():
  return "</html>"

def PrettyHTMLTableHeader(headings):
  string = """
  <H4>Credit for classes / examples is distributed equally among the authors of the class / example</H4>
  <H4>All columns of the table can be sorted by clicking on the appropriate heading</H4>
  <table tableborder="0" cellpadding="1" cellspacing="3" class="sortable">
    <thead>
      <tr>
"""
  for header in headings:
    string += '        <th scope="col" bgcolor="#D5CCB1"><font color="#3300FF" style="bold">' + header + '</th>\n'
  string += """      </tr>
    </thead>
    
    <tbody>"""
  return string

def HTMLTableFooter():
  return '  </tbody>\n</table>\n'

row_colors = ["FAF0D4", "F3DFA8"]
def PrettyHTMLTableRow(row, counter):
  string = '      <tr bgcolor="#' + row_colors[ counter % 2] + '">\n'
  for item in row:
    string += '        <td>' + item + '</td>\n'
  string += '      </tr>'
  return string

def HtmlTableFromTwikiTable(heading, subheadings, header_row, lines):
  if len(header_row) == 0:
    return ""
  rows = ["<h2>" + heading + "</h2>"]
  rows.extend([ "<h3>" + subheading + "</h3>" for subheading in subheadings])
  rows.append(PrettyHTMLTableHeader([ item.strip(' *') for item in header_row.split('|') if len(item) > 0]))
  rows.extend([ PrettyHTMLTableRow([ item.strip(' ') for item in lines[i].split('|') if len(item) > 0], i) for i in xrange(len(lines)) if len(lines[i]) > 0])
  rows.append(HTMLTableFooter())
  return '\n'.join(rows)
