import sys
from collections import deque

def build_graph(sequences, alphabet):
    """
    Phase 1: Efficient Graph Construction.
    Builds an adjacency list where edges connect sequences 
    with a Hamming distance of exactly 1.
    """
    n = len(sequences)
    if n == 0:
        return []
    
    m = len(sequences[0])
    
    # 1. Data structures for fast lookups
    # Set for O(1) existence check
    string_set = set(sequences)
    
    # Dictionary to map string to its index O(1)
    string_to_index = {seq: i for i, seq in enumerate(sequences)}
    
    # 2. Adjacency list
    graph = [[] for _ in range(n)]
    
    # 3. Build edges
    # Instead of O(N^2) comparison, we generate all possible neighbors
    # for each sequence and check if they exist in our dataset.
    for i, seq in enumerate(sequences):
        original_seq_list = list(seq)
        
        for j in range(m): # Iterate over each position
            original_char = original_seq_list[j]
            
            for char in alphabet: # Iterate over alphabet
                if char == original_char:
                    continue
                
                # Create a potential neighbor (Hamming distance = 1)
                original_seq_list[j] = char
                neighbor_str = "".join(original_seq_list)
                
                # Check if this neighbor exists in the input data
                if neighbor_str in string_set:
                    neighbor_index = string_to_index[neighbor_str]
                    
                    # Add edge. We add it here, and when the loop eventually 
                    # reaches the neighbor_index, it will add the back-link to 'i'.
                    graph[i].append(neighbor_index)

            # Restore the character for the next iteration
            original_seq_list[j] = original_char
            
    return graph

def find_connected_components(n, graph):
    """
    Phase 2: Find connected components using BFS.
    """
    visited = [False] * n
    all_components = []

    for i in range(n):
        if not visited[i]:
            # Found a node not yet visited -> start of a new component
            current_component = []
            q = deque([i]) 
            visited[i] = True
            
            while q:
                current_node = q.popleft()
                current_component.append(current_node)
                
                # Traverse neighbors
                for neighbor in graph[current_node]:
                    if not visited[neighbor]:
                        visited[neighbor] = True
                        q.append(neighbor)
            
            # Optional: Sort indices for consistent output
            current_component.sort()
            all_components.append(current_component)
            
    return all_components

def analyze_sequences(sequences, alphabet):
    """
    Main analysis logic and output formatting.
    """
    if not sequences:
        print("Sequence list is empty.")
        return

    print(f"Processing {len(sequences)} sequences...")
    
    # Phase 1: Build Graph
    graph = build_graph(sequences, alphabet)
    
    # Phase 2: Find Components
    components = find_connected_components(len(sequences), graph)
    
    # Output Results
    
    # 1) Number of independent connected components
    print(f"1) Number of connected components: {len(components)}")
    
    # 2) Number of sequences in each component
    # We sort the list of sizes for better readability, 
    # or you can keep them in order of discovery.
    component_sizes = [len(c) for c in components]
    print(f"2) Component sizes: {component_sizes}")
    
    # 3) List of indices for each component
    print(f"3) Indices per component:")
    for i, component in enumerate(components):
        print(f"   Component {i+1}: {component}")

# --- Main Execution Block ---

ALPHABET = "ACDEFGHIKLMNPQRSTVWY*"

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python find_max_subset.py <landscape_file>")
        sys.exit(1)
        
    fname = sys.argv[1]
    genotypes = []
    
    try:
        with open(fname, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # Take the first token as the sequence
                genotypes.append(line.split()[0])
                
        analyze_sequences(genotypes, ALPHABET)
        
    except FileNotFoundError:
        print(f"Error: File '{fname}' not found.")
        sys.exit(1)
