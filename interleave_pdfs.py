from PyPDF2 import PdfReader, PdfWriter

# List your PDF file paths here
pdf_files = ["file1.pdf", "file2.pdf", "file3.pdf", "file4.pdf", "file5.pdf", "file6.pdf"]

# Load all PDFs
readers = [PdfReader(open(pdf, "rb")) for pdf in pdf_files]

# Get the max number of pages across all PDFs
max_pages = max(len(reader.pages) for reader in readers)

writer = PdfWriter()

# Interleave pages
for i in range(max_pages):
    for reader in readers:
        if i < len(reader.pages):
            writer.add_page(reader.pages[i])

# Write to output file
with open("interleaved_output.pdf", "wb") as out_file:
    writer.write(out_file)

print("Interleaved PDF created: interleaved_output.pdf")