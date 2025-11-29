#!/usr/bin/env python3
"""
Test script to verify OpenAI API key is working.
"""

import os
import sys

def test_openai_key():
    """Test if OpenAI API key is set and working."""
    
    # Check if key is set
    api_key = os.getenv("OPENAI_API_KEY")
    
    if not api_key:
        print("❌ ERROR: OPENAI_API_KEY environment variable is NOT set")
        print("\nTo set it, run:")
        print("  export OPENAI_API_KEY='sk-your-api-key-here'")
        print("\nOr add to ~/.bashrc for permanent setup:")
        print("  echo 'export OPENAI_API_KEY=\"sk-your-api-key-here\"' >> ~/.bashrc")
        print("  source ~/.bashrc")
        return False
    
    print(f"✓ OPENAI_API_KEY is set (length: {len(api_key)} characters)")
    print(f"  Key starts with: {api_key[:7]}...")
    
    # Check key format
    if not api_key.startswith("sk-"):
        print("⚠️  WARNING: API key doesn't start with 'sk-' (might be invalid format)")
    
    # Try to import OpenAI
    try:
        from openai import OpenAI
        print("✓ OpenAI library is installed")
    except ImportError:
        print("❌ ERROR: OpenAI library is not installed")
        print("  Install with: pip install openai")
        return False
    
    # Test API connection
    print("\nTesting API connection...")
    try:
        client = OpenAI(api_key=api_key)
        
        # Make a simple test call
        response = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {"role": "system", "content": "You are a helpful assistant."},
                {"role": "user", "content": "Say 'API key is working' if you can read this."}
            ],
            max_tokens=10
        )
        
        result = response.choices[0].message.content
        print(f"✓ API call successful!")
        print(f"  Response: {result}")
        print(f"  Model used: {response.model}")
        print(f"  Tokens used: {response.usage.total_tokens}")
        
        print("\n" + "="*60)
        print("✅ SUCCESS: OpenAI API key is working correctly!")
        print("="*60)
        return True
        
    except Exception as e:
        error_msg = str(e)
        print(f"❌ ERROR: API call failed")
        print(f"  Error: {error_msg}")
        
        if "Invalid API key" in error_msg or "Incorrect API key" in error_msg:
            print("\n  → Your API key appears to be invalid")
            print("  → Check that you copied the full key correctly")
            print("  → Get a new key from: https://platform.openai.com/api-keys")
        elif "You exceeded your current quota" in error_msg or "insufficient_quota" in error_msg:
            print("\n  → Your OpenAI account has insufficient credits/quota")
            print("  → Add credits at: https://platform.openai.com/account/billing")
        elif "Rate limit" in error_msg.lower():
            print("\n  → You've hit the rate limit")
            print("  → Wait a moment and try again")
        else:
            print("\n  → Check your internet connection")
            print("  → Verify your OpenAI account is active")
        
        return False

if __name__ == "__main__":
    success = test_openai_key()
    sys.exit(0 if success else 1)

