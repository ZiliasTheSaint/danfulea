package danfulea.utils;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.pdmodel.PDPageContentStream;
import org.apache.pdfbox.pdmodel.PDPageTree;
import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.apache.pdfbox.pdmodel.font.PDFont;
import org.apache.pdfbox.pdmodel.font.PDType1Font;

import danfulea.phys.egs.EGS4;

/**
 * Utility class for creating PDF file, mainly used for converting long string
 * (which also populate, say, a JTextArea in main application) into a PDF.
 * 
 * @author Dan Fulea, 10 AUG. 2016
 *
 */
@SuppressWarnings("unused")
public class PDFRenderer implements AutoCloseable {
	private final PDDocument doc;

	private PDPageContentStream content = null;
	private float textRenderingLineY = 0;

	/**
	 * Constructor, initializes the PDF document
	 * 
	 * @param doc
	 *            the document
	 */
	public PDFRenderer(PDDocument doc) {
		this.doc = doc;
	}

	/**
	 * Fill the PDF with the text supplied as argument. It uses by default
	 * Helvetica, size 12 font and all margins taken into account are of 1 cm.
	 * If text size is too large then it will be displayed only a part of it
	 * which fits in a row. It does not create new page! 
	 * <p>
	 * This is only useful to print small texts. USE THIS ONLY FOR TESTING.
	 * </p>
	 * 
	 * @param Info
	 *            the String
	 * @throws IOException
	 *             can throw this exception
	 */
	public void renderText(String Info) throws IOException {
		float fontSize = 12f;
		float leading = 1.5f * fontSize;

		float margin = 72;// 1 inch
		margin = (float) (margin / 2.54);// 1 cm

		// footer------------------------
		float footerSize = 10;
		float newPageCondition = leading + margin + 1.5f * footerSize;
		// ------------------------------
		if (content == null || textRenderingLineY < newPageCondition)
			newPage();
		
		textRenderingLineY -= leading;
		
		PDFont fontPlain = PDType1Font.HELVETICA;// .TIMES_ROMAN;//.HELVETICA;//

		PDPage page = doc.getPage(0);// first page
		PDRectangle mediabox = page.getMediaBox();// A4
		float startX = mediabox.getLowerLeftX() + margin;

		content.beginText();
		content.setFont(fontPlain, fontSize);
		
		content.newLineAtOffset(startX, textRenderingLineY);;
		content.showText(Info);
		content.endText();
	}

	/**
	 * Fill the PDF with the text supplied as argument. It uses by default
	 * Helvetica, size 12 font and all margins taken into account are of 1 cm.
	 * The text is splitted into lines based on the presence of the new line
	 * character. It does not create a new page!!
	 * <p>
	 * This is useful to print small texts having new-line characters ("\n") such
	 * those which also populate GUI components (e.g. JTextArea). USE THIS ONLY FOR TESTING.
	 * </p>
	 * 
	 * @param Info
	 *            the String
	 * @throws IOException
	 *             can throw this exception
	 */
	public void renderTextHavingNewLine(String Info) throws IOException {
		float fontSize = 12f;
		float leading = 1.5f * fontSize;

		float margin = 72;// 1 inch
		margin = (float) (margin / 2.54);// 1 cm

		// footer------------------------
		float footerSize = 10;
		float newPageCondition = leading + margin + 1.5f * footerSize;
		// ------------------------------
		if (content == null || textRenderingLineY < newPageCondition)
			newPage();

		textRenderingLineY -= leading;

		PDFont fontPlain = PDType1Font.HELVETICA;// .TIMES_ROMAN;//.HELVETICA;//

		PDPage page = doc.getPage(0);// first page
		PDRectangle mediabox = page.getMediaBox();// A4
		float startX = mediabox.getLowerLeftX() + margin;

		List<String> lines = getLinesFromString(Info);
		
		content.beginText();
		content.setFont(fontPlain, fontSize);
		content.newLineAtOffset(startX, textRenderingLineY);
		for (String line : lines) {
			
			content.showText(line);
			if (!line.equals(lines.get(lines.size() - 1)))// not last line
				textRenderingLineY -= leading;

			content.newLineAtOffset(0, -leading);
			// Move to the start of the next line, offset from the start of
			// the CURRENT line by (tx, ty).
		}
		
		content.endText();
	}

	/**
	 * Sets the title on top of the page center aligned. If title is too large,
	 * then it is displayed on more lines not centered. It uses helvetica bold
	 * font of default size 12 and all margins taken into account are of 1 cm.
	 * 
	 * @param title
	 *            the title
	 * @throws IOException
	 *             can throw this exception
	 */
	public void setTitle(String title) throws IOException {
		
		if (content == null)
			newPage();

		float fontSize = titleFontSize;//16f;
		float leading = 1.5f * fontSize;

		PDFont font = PDType1Font.HELVETICA_BOLD;// .TIMES_ROMAN;//.HELVETICA;
		float titleWidth = font.getStringWidth(title) / 1000f * fontSize;
		
		PDPage page = doc.getPage(0);// first page
		PDRectangle mediabox = page.getMediaBox();// A4
		
		// From the rectangle you can grab the width and height and convert it
		// to
		// your purpose. The 72 is the dpi value that you divide and optional
		// you
		// can multiply it by 2.54 to obtain cm. so getWidth=value in PostScript
		// points
		// To obtain inch =>value/72. To obtain cm =>value/72*2.54

		float margin = 72;// 1 inch
		margin = (float) (margin / 2.54);// 1 cm
		float width = mediabox.getWidth() - 2 * margin;

		float xP = margin + (width - titleWidth) / 2;// this centers on screen!!!

		float startX = mediabox.getLowerLeftX() + margin;
		
		// in case text is too long
		List<String> lines = getLines(title, fontSize, font, width);
		// end text too long=============

		content.beginText();
		content.setFont(font, fontSize);
		
		textRenderingLineY -= leading;//some space
		
		if (lines.size() > 1)
			content.newLineAtOffset(startX, textRenderingLineY);
		else
			content.newLineAtOffset(xP, textRenderingLineY);

		for (String line : lines) {

			content.showText(line);
			if (!line.equals(lines.get(lines.size() - 1)))// not last line
				textRenderingLineY -= leading;

			if (lines.size() > 1) {
				
				content.newLineAtOffset(0, -leading);
				// Move to the start of the next line, offset from the start of
				// the CURRENT line by (tx, ty).
			} else
				content.newLineAtOffset(xP, textRenderingLineY);
		}
		
		content.endText();

	}

	/**
	 * Sets the subtitle on top of the page below the title center aligned. If
	 * subtitle is too large, then it is displayed on more lines not centered.
	 * It uses helvetica bold font of default size 10 and all margins taken into account
	 * are of 1 cm.
	 * 
	 * @param subtitle
	 *            the subtitle
	 * @throws IOException
	 *             can throw this exception
	 */
	public void setSubTitle(String subtitle) throws IOException {
		
		if (content == null)
			newPage();

		float fontSize = subtitleFontSize;//14f;
		float leading = 1.5f * fontSize;

		textRenderingLineY -= leading;//some space
		textRenderingLineY -= leading;

		PDFont font = PDType1Font.HELVETICA_BOLD;// .TIMES_ROMAN;//.HELVETICA;
		float titleWidth = font.getStringWidth(subtitle) / 1000f * fontSize;

		PDPage page = doc.getPage(0);// first page
		PDRectangle mediabox = page.getMediaBox();// A4

		float margin = 72;// 1 inch
		margin = (float) (margin / 2.54);// 1 cm
		float width = mediabox.getWidth() - 2 * margin;

		float xP = margin + (width - titleWidth) / 2;
		float startX = mediabox.getLowerLeftX() + margin;

		// in case text is too long
		List<String> lines = getLines(subtitle, fontSize, font, width);
		// end text too long=============

		PDFont fontPlain = PDType1Font.HELVETICA_BOLD;// .TIMES_ROMAN;//.HELVETICA;
		content.beginText();
		content.setFont(fontPlain, fontSize);

		if (lines.size() > 1)
			content.newLineAtOffset(startX, textRenderingLineY);
		else
			content.newLineAtOffset(xP, textRenderingLineY);
		
		for (String line : lines) {
			
			content.showText(line);
			if (!line.equals(lines.get(lines.size() - 1)))// not last line
				textRenderingLineY -= leading;

			if (lines.size() > 1) {
				
				content.newLineAtOffset(0, -leading);
				// Move to the start of the next line, offset from the start of
				// the CURRENT line by (tx, ty).
			} else
				content.newLineAtOffset(xP, textRenderingLineY);
		}
		
		content.endText();

	}

	/**
	 * Get the lines arrays based on input text. If text doesn't fit the row,
	 * then it is split into lines using space (" ") character as splitter.
	 * 
	 * @param text
	 *            the text
	 * @param fontSize
	 *            text font size
	 * @param pdfFont
	 *            text font type
	 * @param width
	 *            the width limit to fit into
	 * @return the line array
	 * @throws IOException
	 *             can throw this exception
	 */
	private List<String> getLines(String text, float fontSize, PDFont pdfFont, float width) throws IOException {
		List<String> lines = new ArrayList<String>();

		int lastSpace = -1;
		while (text.length() > 0) {
			int spaceIndex = text.indexOf(' ', lastSpace + 1);
			// Returns the index within this string of the first occurrence
			// of the specified substring, starting at the specified index.
			if (spaceIndex < 0)
				spaceIndex = text.length();
			String subString = text.substring(0, spaceIndex);
			float size = fontSize * pdfFont.getStringWidth(subString) / 1000f;
	
			if (size > width) {
				if (lastSpace < 0)
					lastSpace = spaceIndex;
				subString = text.substring(0, lastSpace);
				lines.add(subString);
				text = text.substring(lastSpace).trim();
				// the substring begins with the character at the specified
				// index and extends to the end of this string.
		
				lastSpace = -1;
			} else if (spaceIndex == text.length()) {
				lines.add(text);
				text = "";
			} else {
				lastSpace = spaceIndex;
			}
		}

		return lines;
	}

	/**
	 * Get the lines arrays based on input text. The text is split into lines
	 * using new line ("\n") character as splitter.
	 * 
	 * @param text
	 *            the text
	 * @return the line array
	 * @throws IOException
	 *             can throw this exception
	 */
	private List<String> getLinesFromString(String text) throws IOException {
		List<String> lines = new ArrayList<String>();

		int lastSpace = -1;
		while (text.length() > 0) {
			int spaceIndex = text.indexOf('\n', lastSpace + 1);
			// Returns the index within this string of the first occurrence
			// of the specified substring, starting at the specified index.
			if (spaceIndex < 0)
				spaceIndex = text.length();
			String subString = text.substring(0, spaceIndex);
	
			// make split condition
			if (spaceIndex != text.length()) {
				if (lastSpace < 0)
					lastSpace = spaceIndex;
				subString = text.substring(0, lastSpace);
				lines.add(subString);
				text = text.substring(lastSpace).trim();
				// the substring begins with the character at the specified
				// index and extends to the end of this string.
			
				lastSpace = -1;
			} else if (spaceIndex == text.length()) {
				lines.add(text);
				text = "";
			} 
		}

		return lines;
	}

	/**
	 * Creates new A4 page.
	 * 
	 * @throws IOException
	 *             can throw this exception
	 */
	private void newPage() throws IOException {
		close();

		PDPage page = new PDPage(PDRectangle.A4);// A4
		doc.addPage(page);
		content = new PDPageContentStream(doc, page);
		content.setNonStrokingColor(Color.BLACK);

		PDRectangle mediabox = page.getMediaBox();// A4
		float margin = 72;// 1 inch
		margin = (float) (margin / 2.54);// 1 cm

		textRenderingLineY = mediabox.getUpperRightY() - margin;
	}

	/**
	 * Closes the current page
	 */	
	public void close() throws IOException {
		if (content != null) {
			content.close();
			content = null;
		}
	}

	/**
	 * Sets the page number on all pages center aligned. It uses helvetica bold
	 * font of size equal with regularFontSize and all margins taken into account are of 1 cm.
	 * 
	 * @throws IOException
	 *             can throw this exception
	 */
	public void addPageNumber() throws IOException {
	
		PDPageTree allPages = doc.getPages();
		PDFont font = PDType1Font.HELVETICA_BOLD;
		float fontSize = regularFontSize;//10.0f;

		float margin = 72;// 1 inch
		margin = (float) (margin / 2.54);// 1 cm

		for (int i = 0; i < allPages.getCount(); i++) {
			PDPage page = (PDPage) allPages.get(i);
			PDRectangle pageSize = page.getMediaBox();

			// -----------------------------------
			int pageNumber = i + 1;
			String message = "" + pageNumber;
			// ----------------------------------
			float stringWidth = font.getStringWidth(message) * fontSize / 1000f;

			// --------------------------
			float startY = pageSize.getLowerLeftY() + margin;
			float width = pageSize.getWidth() - 2 * margin;
			float startX = margin + (width - stringWidth) / 2;
			// -----------------------------

			// calculate to center of the page
			// int rotation = page.getRotation(); // page.getr
			// boolean rotate = rotation == 90 || rotation == 270;
			// float pageWidth = rotate ? pageSize.getHeight() :
			// pageSize.getWidth();
			// float pageHeight = rotate ? pageSize.getWidth() :
			// pageSize.getHeight();

			// float centeredXPosition = rotate ? pageHeight / 2f : (pageWidth -
			// stringWidth) / 2f;
			// float centeredYPosition = rotate ? (pageWidth - stringWidth) / 2f
			// : pageHeight / 2f;
			// append the content to the existing stream
			PDPageContentStream contentStream = // new PDPageContentStream(doc,
												// page, true, true,true);
					new PDPageContentStream(doc, page, PDPageContentStream.AppendMode.APPEND, true, true);
			contentStream.beginText();
			// set font and font size
			contentStream.setFont(font, fontSize);
			// set text color to red
			contentStream.setNonStrokingColor(Color.BLACK);// 255, 0, 0);
			// ---------------
			contentStream.newLineAtOffset(startX, startY);// at footer
			/// --------------
			/*
			 * if (rotate) { // rotate the text according to the page rotation
			 * // contentStream.setTextRotation(Math.PI/2, centeredXPosition, //
			 * centeredYPosition); contentStream
			 * .setTextMatrix(Matrix.getRotateInstance(Math.PI / 2,
			 * centeredXPosition, centeredYPosition)); } else {
			 * //contentStream.setTextTranslation(centeredXPosition,
			 * centeredYPosition);
			 * contentStream.setTextMatrix(Matrix.getTranslateInstance(
			 * centeredXPosition, centeredYPosition)); }
			 */
			contentStream.showText(message);// drawString(message);
			contentStream.endText();
			contentStream.close();
		}

		// doc.save( outfile );
		// }
		// finally
		// {
		// if( doc != null )
		// {
		// doc.close();
		// }
		// }
		
	}
	
	private float subtitleFontSize = 10f;
	/**
	 * Setup subtitle font size. Default is 10f; f comes from float!
	 * @param subtitleFontSize subtitleFontSize
	 */
	public void setSubtitleFontSize(float subtitleFontSize){
		this.subtitleFontSize=subtitleFontSize;
	}
	
	private float titleFontSize = 12f;
	/**
	 * Setup title font size. Default is 12f; f comes from float!
	 * @param titleFontSize titleFontSize
	 */
	public void setTitleFontSize(float titleFontSize){
		this.titleFontSize=titleFontSize;
	}
	
	private float regularFontSize = 7f;
	/**
	 * Setup regular font size. Default is 7f; f comes from float!
	 * @param regularFontSize regularFontSize
	 */
	public void setRegularFontSize(float regularFontSize){
		this.regularFontSize=regularFontSize;
	}
	
	/** Fill the PDF with the text supplied as argument. It uses by default
	 * COURIER, default font size of 7 (which can be changed via setRegularFontSize) and all margins taken into account are of 1 cm. 
	 * The text is splitted into lines based on the presence of the new line
	 * character. It does create a new page if needed! Also if a line doesn't fit 
	 * the page width, a new line is automatically created from last SPACE character.
	 * <p>
	 * This is useful to print all kind of texts, for instance 
	 * those texts which also populate GUI components (e.g. JTextArea).
	 * <p>
	 * Note that the text must be trimmed and TAB-free in order to be successfully displayed (this is handled automatically by this routine). Otherwise, 
	 * this version of PDFBox library complains about u0009 invalid font type when encounter TAB and 
	 * TAB-like characters. This may affects texts which uses those characters, not a really big issue since they are still displayed reasonably well, 
	 * but it is a known shortcoming for now. iText library performs better in this regard but it is GPL-ed, so I try to avoid any restrictive 
	 * softwares (i.e. force to create only free software!). 
	 * 
	 * @param text
	 *            the String
	 * @throws IOException
	 *             can throw this exception
	 */
	public void renderTextEnhanced(String text) throws IOException{
		float fontSize = regularFontSize;//14f;
		float leading = 1.5f * fontSize;

		float margin = 72;// 1 inch
		margin = (float) (margin / 2.54);// 1 cm
		// footer------------------------
		float footerSize = 10;
		float newPageCondition = leading + margin + 1.5f * footerSize;
		// ------------------------------
		//PDType1Font.COURIER
		PDFont fontPlain = PDType1Font.COURIER;// .TIMES_ROMAN;//.HELVETICA;//
		PDPage page = doc.getPage(0);// first page
		PDRectangle mediabox = page.getMediaBox();// A4
		float startX = mediabox.getLowerLeftX() + margin;
		float width = mediabox.getWidth() - 2 * margin;
		
		
		if (content == null || textRenderingLineY < newPageCondition)
			newPage();

		textRenderingLineY -= leading;//some space
		textRenderingLineY -= leading;

		content.beginText();
		content.setFont(fontPlain, fontSize);
		content.newLineAtOffset(startX, textRenderingLineY);		
		//=========
		String[] wholelines = text.split("\n");//looking for EOL character
		
		for (int i = 0; i<wholelines.length;i++){//37;i++){//wholelines.length;i++){
			//-------------------------
			wholelines[i] = wholelines[i].trim();//trim is vital for removing tabs...
			//if however wholelines[i] contains TABS inside they are not removed by trim. 
			//ALso, replace and replaceAll do not WORK!=>use split.
			String[] splitted = wholelines[i].split("\t");//;looking for TABS!!!
			String finalStr="";
			for (int k = 0;k<splitted.length;k++){
				finalStr = finalStr+" "+splitted[k];
			}
			//finalStr = wholelines[i];
			//---------------
			//List<String> linesOfLine = getLines(wholelines[i], fontSize, fontPlain, width);
			List<String> linesOfLine = getLines(finalStr, fontSize, fontPlain, width);
			for (int j=0;j<linesOfLine.size();j++){
				//lines.add(linesOfLine.get(j));
				//renderLine(linesOfLine.get(j), newPageCondition, leading, startX);
				if (content == null || textRenderingLineY < newPageCondition){
					newPage();//also set textRenderingLineY
		
					content.beginText();
					content.setFont(fontPlain, fontSize);
					content.newLineAtOffset(startX, textRenderingLineY);
				}
				
				textRenderingLineY -= leading;
	
				String strng = linesOfLine.get(j);
								
				content.showText(strng);//finalStr);//strng);//"test");//wholelines[i]);//
				content.newLineAtOffset(0, -leading);	
			}
		}
		
		//================
		content.endText();
	}	
	// =================================================TESTING
		/*public static void main(String[] args) {
			PDDocument doc = new PDDocument();
			PDFRenderer renderer = new PDFRenderer(doc);		
			
			try{
				renderer.setTitle("MyApp really  very very very very very very long title which is generated Report");//, marginTop);
				renderer.setSubTitle("Report generated on May 12, 2020");//, marginTop);
				for (int i = 0; i < 2000; i++)
				{
			    
					//renderer.renderText("hello" + i);//, margin);//60);
				
				}//WORKS!!!!!!!!
				
				String str = "This is "+"\n";
				str= str+ "abc"+"\n";
				str= str+ "The most powerful app"+"\n";
				str= str+ "on earth";
			
				str ="		   "
							+ "mu_tr/rho = <E*mu_tr/rho>/Eave     = "
						+ EGS4.format(12.45, 12, false)+ " 	cm^2/g " + "	+/- "
							+ EGS4.format(100 * 1.2547, 12) + " %";;
				//String[] splitted = str.split("\t");//;looking for TABS!!!
				//String finalStr="";
				//for (int k = 0;k<splitted.length;k++){
				//	finalStr = finalStr+" "+splitted[k];
				//}
				//renderer.renderTextHavingNewLine(str);//works!!!!!!!!!!!!!!!!
				renderer.renderTextEnhanced(str);//finalStr);
				
				renderer.addPageNumber();
				renderer.close();		
				doc.save(new File("renderSimple.pdf"));
				doc.close();
			}catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}*/
		
		
	// ================================================END TESTING
}
