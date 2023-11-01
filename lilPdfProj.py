import pyttsx3
import PyPDF2

def read_pdf_with_speed(file_path, speed):
    # Open the PDF file
    pdf_file = open(file_path, 'rb')

    # Create a PdfFileReader object
    pdf_reader = PyPDF2.PdfFileReader(pdf_file)

    # Initialize the text-to-speech engine
    engine = pyttsx3.init()

    # Set the speech speed
    rate = engine.getProperty('rate')
    engine.setProperty('rate', rate + speed)

    # Read each page and convert the text to speech
    for page_num in range(pdf_reader.numPages):
        page = pdf_reader.getPage(page_num)
        text = page.extractText()
        engine.say(text)
        engine.runAndWait()

    # Close the PDF file
    pdf_file.close()

# Use the function
read_pdf_with_speed('your_pdf_file.pdf', 50)
