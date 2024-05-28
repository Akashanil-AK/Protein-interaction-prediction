import streamlit as st
import joblib
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np

# Load your trained model
model = joblib.load ('Final Model.pkl')
top_n_features = joblib.load('top_n_features.pkl')

def calculate_aac(sequence):
    protein_analysis = ProteinAnalysis(sequence)
    aac = protein_analysis.get_amino_acids_percent()
    return aac

def calculate_k_spaced_pairs(sequence, k):
    pairs = {}
    length = len(sequence)
    for i in range(length - k):
        pair = sequence[i] + sequence[i + k]
        pairs[pair] = pairs.get(pair, 0) + 1
    return pairs

def calculate_pse_aac_type_ii(sequence, k, max_length):
    aac = calculate_aac(sequence)
    pairs = calculate_k_spaced_pairs(sequence, k)
    dipeptides = [sequence[i:i+2] for i in range(len(sequence)-1)]
    dipeptide_counts = {dipeptide: dipeptides.count(dipeptide) for dipeptide in set(dipeptides)}
    pse_aac = list(aac.values())
    pse_aac.extend(pairs.values())
    pse_aac.extend(dipeptide_counts.values())
    if len(pse_aac) < max_length:
        pse_aac.extend([0] * (max_length - len(pse_aac)))
    elif len(pse_aac) > max_length:
        pse_aac = pse_aac[:max_length]
    return np.array(pse_aac)

# Preprocess a single input sequence
def preprocess_sequence(sequence1, sequence2, k=3, max_length=3601):
    pse_aac_1 = calculate_pse_aac_type_ii(sequence1, k, max_length)
    pse_aac_2 = calculate_pse_aac_type_ii(sequence2, k, max_length)
    combined_features = np.concatenate((pse_aac_1, pse_aac_2))
    return combined_features

# Predict function for a single input
def predict_interaction(sequence1, sequence2):
    # Preprocess the input sequences
    input_features = preprocess_sequence(sequence1, sequence2)
    input_features = input_features.reshape(1, -1)  # Reshape for a single sample

    # Select top N features
    input_features_selected = input_features[:, top_n_features]

st.title('Protein Interaction Predictor')

st.write("""
Input the protein features to predict interactions.
""")

# Input fields for protein features
sequence1 = st.text_input('Feature 1')
sequence2 = st.text_input('Feature 2')
# Add more features as needed



if st.button('Predict'):
    features = [str (sequenc1), str (sequence2)]  # Extend this list with more features
    input_data = preprocess_sequence(features)
    input_select = top_n_features(input_data)
    prediction = model.predict(input_select)
    if prediction[0] == 1:
            st.write('The given sequences interact.')
    else:
            st.write('No interaction between the given sequences.')
    
    
