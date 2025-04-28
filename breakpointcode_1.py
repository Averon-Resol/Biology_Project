import pandas as pd
import os
import random
import networkx as nx
import matplotlib.pyplot as plt

# =========================
# Step 1: Read GTF and Build Genome1
# =========================

# Get current working folder
current_folder = os.getcwd()

# GTF filename (keep GTF file and script in same folder)
gtf_filename = "gencode.v47.basic.annotation.gtf"   # Change if needed
gtf_file = os.path.join(current_folder, gtf_filename)

genes = []

# Read GTF file line by line
with open(gtf_file, 'r') as file:
    for line in file:
        if line.startswith("#"):
            continue
        fields = line.strip().split('\t')
        
        if len(fields) < 9:
            continue
        
        feature_type = fields[2]
        
        if feature_type == 'gene':
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            
            gene_id = None
            for entry in attributes.split(";"):
                if "gene_id" in entry:
                    gene_id = entry.strip().split(' ')[1].replace('"','')
                    break
            
            if gene_id:
                genes.append((chrom, start, strand, gene_id))

# Build DataFrame
df = pd.DataFrame(genes, columns=['chrom', 'start', 'strand', 'gene_id'])

# Optional: filter to 1 chromosome for small test
# df = df[df['chrom'] == 'chr1']

# Sort by chromosome and start
df = df.sort_values(by=['chrom', 'start'])

# Assign number to each gene
gene_to_num = {gene: idx+1 for idx, gene in enumerate(df['gene_id'])}

# Build signed genome1
genome1 = []
for idx, row in df.iterrows():
    gene_num = gene_to_num[row['gene_id']]
    if row['strand'] == '+':
        genome1.append(gene_num)
    else:
        genome1.append(-gene_num)

print("\n✅ Signed Genome1 (first 100 genes):")
print(genome1[:100])

# Save Genome1 to file
with open('signed_genome1_list.txt', 'w') as f:
    f.write(' '.join(map(str, genome1)))

print("\nSigned Genome1 saved as 'signed_genome1_list.txt'.")

# =========================
# Step 2: Simulate Mutation to Create Genome2
# =========================

# Copy genome1
genome2 = genome1.copy()

# Pick a random region to invert
start = random.randint(0, len(genome2) - 20)  # Start point
end = start + random.randint(5, 10)            # Length of 5 to 10 genes

# Invert the region
genome2[start:end] = [-x for x in reversed(genome2[start:end])]

print(f"\n✅ Simulated Mutation: Inverted genes from index {start} to {end-1}.")

print("\n✅ Signed Genome2 (first 100 genes):")
print(genome2[:100])

# Save Genome2 to file
with open('signed_genome2_list.txt', 'w') as f:
    f.write(' '.join(map(str, genome2)))

print("\nSigned Genome2 saved as 'signed_genome2_list.txt'.")

# =========================
# Step 3: Compare Genomes and Count Breakpoints
# =========================

def count_breakpoints(genome1, genome2):
    """
    Count the number of breakpoints between two genomes
    based on adjacent gene pairs.
    """
    # Build adjacency sets
    adj1 = set((genome1[i], genome1[i+1]) for i in range(len(genome1)-1))
    adj2 = set((genome2[i], genome2[i+1]) for i in range(len(genome2)-1))
    
    # Breakpoints = places where adjacencies differ
    breakpoints = adj1.symmetric_difference(adj2)
    
    return len(breakpoints)

def find_breakpoints(genome1, genome2):
    """
    Find and list the breakpoints between two genomes
    based on adjacency differences.
    """
    adj1 = set((genome1[i], genome1[i+1]) for i in range(len(genome1)-1))
    adj2 = set((genome2[i], genome2[i+1]) for i in range(len(genome2)-1))
    
    breakpoints = adj1.symmetric_difference(adj2)
    
    print(f"\n✅ Total number of Breakpoints: {len(breakpoints)}\n")
    print("Breakpoints (gene pairs that are broken):")
    for b in breakpoints:
        print(b)

# Call to find and display breakpoints
find_breakpoints(genome1, genome2)

# Just before draw
small_genome1 = genome1[:200]
small_genome2 = genome2[:200]


def draw_breakpoint_graph(small_genome1, small_genome2):
    """
    Draw the breakpoint graph between two genomes with clearer distinction.
    """
    G = nx.Graph()
    
    # Create directed adjacencies (order matters)
    adj1 = [(small_genome1[i], small_genome1[i+1]) for i in range(len(small_genome1)-1)]
    adj2 = [(small_genome2[i], small_genome2[i+1]) for i in range(len(small_genome2)-1)]
    
    # Convert to hashable frozensets for comparison (ignoring direction)
    adj1_set = {frozenset([u, v]) for u, v in adj1}
    adj2_set = {frozenset([u, v]) for u, v in adj2}
    
    # Calculate unique edges
    only_in_genome1 = adj1_set - adj2_set
    only_in_genome2 = adj2_set - adj1_set
    in_both = adj1_set & adj2_set
    
    # Add nodes first (to ensure all nodes are included)
    all_nodes = set()
    for u, v in adj1 + adj2:
        all_nodes.add(u)
        all_nodes.add(v)
    
    for node in all_nodes:
        G.add_node(node)
    
    # Add edges with appropriate colors
    for u, v in adj1:
        edge_set = frozenset([u, v])
        if edge_set in only_in_genome1:
            G.add_edge(u, v, color='red')
    
    for u, v in adj2:
        edge_set = frozenset([u, v])
        if edge_set in only_in_genome2:
            G.add_edge(u, v, color='blue')
    
    for edge_set in in_both:
        u, v = tuple(edge_set)
        G.add_edge(u, v, color='purple')
    
    # Get edge colors for drawing
    colors = [G[u][v]['color'] for u, v in G.edges()]
    
    # Print some debugging info
    print(f"Total nodes: {len(G.nodes())}")
    print(f"Total edges: {len(G.edges())}")
    print(f"Red edges (genome1 only): {colors.count('red')}")
    print(f"Blue edges (genome2 only): {colors.count('blue')}")
    print(f"Purple edges (both): {colors.count('purple')}")
    
    # Use a better layout for this type of graph
    pos = nx.kamada_kawai_layout(G)
    
    plt.figure(figsize=(12, 10))
    nx.draw(G, pos, edge_color=colors, with_labels=True, 
            node_size=300, font_size=8, width=2.0)
    
    # Add a legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='red', lw=2, label='Genome1 Only'),
        Line2D([0], [0], color='blue', lw=2, label='Genome2 Only'),
        Line2D([0], [0], color='purple', lw=2, label='Both Genomes')
    ]
    plt.legend(handles=legend_elements, loc="upper right")
    
    plt.title("Breakpoint Graph")
    plt.savefig('breakpoint_graph.png')
    plt.show()

# Let's also add a function to verify the mutation in genome2
def verify_mutation(genome1, genome2, start_index, end_index):
    """
    Verify that the mutation was correctly applied
    """
    print("\n=== MUTATION VERIFICATION ===")
    print(f"Original segment in genome1: {genome1[start_index:end_index]}")
    print(f"Inverted segment in genome2: {genome2[start_index:end_index]}")
    
    # Check if it's correctly inverted
    original = genome1[start_index:end_index]
    should_be = [-x for x in reversed(original)]
    actual = genome2[start_index:end_index]
    
    if should_be == actual:
        print("✅ Mutation correctly applied!")
    else:
        print("❌ Mutation NOT correctly applied!")
        print(f"Expected: {should_be}")
        print(f"Actual:   {actual}")
    
    # Check for changed adjacencies
    before_mutation = genome1[start_index-1] if start_index > 0 else None
    after_mutation = genome1[end_index] if end_index < len(genome1) else None
    
    print("\nBreakpoints created by mutation:")
    if start_index > 0:
        print(f"Before inversion: {genome1[start_index-1]} -> {genome1[start_index]}")
        print(f"After inversion:  {genome2[start_index-1]} -> {genome2[start_index]}")
    
    if end_index < len(genome1):
        print(f"Before inversion: {genome1[end_index-1]} -> {genome1[end_index]}")
        print(f"After inversion:  {genome2[end_index-1]} -> {genome2[end_index]}")

# To use this verification, add this after your mutation:
verify_mutation(genome1, genome2, start, end)
# Call it
draw_breakpoint_graph(small_genome1, small_genome2)




