# OpenAI Biomarker Extractor Setup Guide

This guide explains how to set up and use the OpenAI-based biomarker extraction feature.

## 📋 Prerequisites

1. **OpenAI API Key**: You need an OpenAI API account and API key
2. **Python Dependencies**: Install required packages
3. **Environment Variable**: Set `OPENAI_API_KEY`

## 🚀 Step-by-Step Setup

### Step 1: Get OpenAI API Key

1. Go to [OpenAI Platform](https://platform.openai.com/)cond
2. Sign up or log in to your account
3. Navigate to **API Keys** section: https://platform.openai.com/api-keys
4. Click **"Create new secret key"**
5. Copy the API key (you'll only see it once!)

### Step 2: Install Required Dependencies

The required packages are already in `requirements.txt`. Install them:

```bash
# Activate your conda environment
conda activate ExoGPT

# Install OpenAI and related packages
pip install openai>=1.0.0 openpyxl>=3.1.0 PyPDF2>=3.0.0 python-docx>=1.1.0
```

Or install all requirements:

```bash
pip install -r requirements.txt
```

### Step 3: Set Environment Variable

You need to set the `OPENAI_API_KEY` environment variable. Choose one method:

#### Option A: Export in Terminal (Temporary - for current session)

```bash
export OPENAI_API_KEY="sk-your-api-key-here"
```

**Note**: This only works for the current terminal session. You'll need to set it again if you open a new terminal.

#### Option B: Add to `.bashrc` or `.zshrc` (Permanent)

```bash
# Add to ~/.bashrc (or ~/.zshrc for zsh)
echo 'export OPENAI_API_KEY="sk-your-api-key-here"' >> ~/.bashrc
source ~/.bashrc
```

#### Option C: Create `.env` File (Recommended for development)

Create a `.env` file in the project root:

```bash
cd /home/david/.cursor-tutor/ExoGPT
echo "OPENAI_API_KEY=sk-your-api-key-here" > .env
```

Then modify `backend/app.py` to load from `.env` (if using python-dotenv):

```bash
pip install python-dotenv
```

Add to `backend/app.py`:
```python
from dotenv import load_dotenv
load_dotenv()
```

#### Option D: Set in System Environment (Linux)

```bash
# Add to /etc/environment (requires sudo)
sudo nano /etc/environment
# Add: OPENAI_API_KEY="sk-your-api-key-here"
```

### Step 4: Verify Setup

Test if the API key is set correctly:

```bash
# Check if environment variable is set
echo $OPENAI_API_KEY

# Or test in Python
python -c "import os; print('API Key set:', bool(os.getenv('OPENAI_API_KEY')))"
```

## 📝 How to Use

### Method 1: Upload via Backend API (Programmatic)

The backend endpoint `/api/step1/publication/upload` will automatically:
1. Save the file to `data/publications/`
2. Extract text from PDF/DOC/TXT
3. Call OpenAI API to extract biomarkers
4. Save results to `data/publications/publications_data.xlsx`

**Example using curl:**
```bash
curl -X POST http://localhost:5000/api/step1/publication/upload \
  -F "file=@/path/to/your/publication.pdf"
```

### Method 2: Manual File Upload (Recommended)

1. **Copy your publication files** directly to the `data/publications/` directory:
   ```bash
   cp /path/to/your/publication.pdf /home/david/.cursor-tutor/ExoGPT/data/publications/
   ```

2. **Run the extraction script** manually:
   ```python
   from exo_gpt.publication_extractor import extract_publication_data
   import os
   
   PROJECT_ROOT = "/home/david/.cursor-tutor/ExoGPT"
   publications_dir = os.path.join(PROJECT_ROOT, "data", "publications")
   file_path = os.path.join(publications_dir, "your_publication.pdf")
   
   # Extract and save to Excel
   result = extract_publication_data(file_path, publications_dir)
   print(f"Extracted {len(result['biomarkers'])} biomarkers")
   ```

3. **Or use the GUI**: The GUI will automatically detect files in `data/publications/` and show them. Click on a publication to see extracted biomarkers.

## 🔍 How It Works

1. **File Detection**: System scans `data/publications/` directory
2. **Excel Check**: Checks if `publications_data.xlsx` already contains data for this file
3. **If Not Found**:
   - Extracts text from PDF/DOC/TXT
   - Sends to OpenAI API (gpt-4o-mini model)
   - Extracts:
     - Publication metadata (title, journal, authors, year, DOI, PMID)
     - Biomarkers (gene symbols, UniProt IDs, logFC, p-values, etc.)
   - Saves to Excel file
4. **If Found**: Uses existing data from Excel (skips API call to save costs)

## 📊 Output Format

Data is saved to `data/publications/publications_data.xlsx` with two sheets:

### Sheet 1: Publications
- filename
- title
- journal
- authors
- year
- doi
- pmid

### Sheet 2: Biomarkers
- publication_filename
- gene_symbol
- uniprot
- analyte_type
- disease
- biofluid
- logfc
- p_value
- fdr
- method

## 💡 Tips

1. **Cost Optimization**: The system checks Excel first, so re-processing the same file won't call OpenAI API again
2. **File Formats**: Supports PDF, DOC, DOCX, and TXT files
3. **Large Files**: Text is truncated to ~100,000 characters (25k tokens) to stay within API limits
4. **View Results**: Use the GUI to click on publications and see extracted biomarkers
5. **Manual Processing**: You can manually trigger extraction by running the Python script

## ⚠️ Troubleshooting

### Error: "OPENAI_API_KEY environment variable not set"
- Make sure you've set the environment variable
- Restart your terminal/IDE after setting it
- Verify with: `echo $OPENAI_API_KEY`

### Error: "Failed to extract text from PDF"
- Check if the PDF is encrypted or corrupted
- Try converting to TXT format first
- Ensure PyPDF2 is installed: `pip install PyPDF2>=3.0.0`

### Error: "OpenAI API error"
- Check your API key is valid
- Verify you have credits in your OpenAI account
- Check API rate limits

### No biomarkers extracted
- The publication might not contain biomarker information
- Try a different publication
- Check the extracted text quality

## 🔐 Security Note

**Never commit your API key to git!**

- Add `.env` to `.gitignore` if using .env file
- Don't share your API key publicly
- Use environment variables instead of hardcoding

## 📚 Additional Resources

- [OpenAI API Documentation](https://platform.openai.com/docs)
- [OpenAI Pricing](https://openai.com/pricing)
- Model used: `gpt-4o-mini` (cost-effective for extraction tasks)

