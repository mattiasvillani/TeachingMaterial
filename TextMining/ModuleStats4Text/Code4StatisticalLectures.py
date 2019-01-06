# Statistical Models for Text

# Import NLTK module
import nltk
# nltk.download() # Uncomment this line if you have installed all the text corpora from NLTK. 
                  # Beware that there will be a new window opened in the background (hidden) where you need to make choices.

# Reading all the books in NLTKs collection to memory (there are nine plus some sentences).
from nltk.book import *

# If you only want one of them, say text1, do the obvious:
from nltk.book import text7

# Define your own text
myText = ["This", "is", "my","text","and","there","is","nothing","you","can","do","about","it","!"]

bigramsText = nltk.bigrams(myText) # This sets up a so called generator of bigrams, ready to produce your bigrams
list(bigramsText)                  # This returns a list with all the bigrams in the text

# POS-tagging
nltk.pos_tag(myText)

nltk.pos_tag(nltk.word_tokenize('I will race you tomorrow')) # Wrongly tags 'race' as a noun (NN)

nltk.pos_tag(nltk.word_tokenize('There is a race tomorrow')) # Rightly tags 'race' as a noun (NN)