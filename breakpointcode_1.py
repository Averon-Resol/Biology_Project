import pandas as pd
import os
import random
import networkx as nx
import matplotlib.pyplot as plt

current_folder = os.getcwd()
gtf_filename = "gencode.v47.basic.annotation.gtf"
gtf_file = os.path.join(current_folder, gtf_filename)

genes = []

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

df = pd.DataFrame(genes, columns=['chrom', 'start', 'strand', 'gene_id'])
df = df.sort_values(by=['chrom', 'start'])
gene_to_num = {gene: idx+1 for idx, gene in enumerate(df['gene_id'])}

genome1 = []
for idx, row in df.iterrows():
    gene_num = gene_to_num[row['gene_id']]
    if row['strand'] == '+':
        genome1.append(gene_num)
    else:
        genome1.append(-gene_num)

print("\n✅ Signed Genome1 (first 100 genes):")
print(genome1[:100])

with open('signed_genome1_list.txt', 'w') as f:
    f.write(' '.join(map(str, genome1)))

print("\nSigned Genome1 saved as 'signed_genome1_list.txt'.")

genome2 = genome1.copy()

# Create a LARGER mutation to make differences more visible
start = random.randint(0, len(genome2) - 50)  # Start point
end = start + random.randint(20, 30)          # Length of 20 to 30 genes

# Invert the region
genome2[start:end] = [-x for x in reversed(genome2[start:end])]

print(f"\n✅ Simulated Mutation: Inverted genes from index {start} to {end-1}.")

print("\n✅ Signed Genome2 (first 100 genes):")
print(genome2[:100])

with open('signed_genome2_list.txt', 'w') as f:
    f.write(' '.join(map(str, genome2)))

print("\nSigned Genome2 saved as 'signed_genome2_list.txt'.")

def count_breakpoints(genome1, genome2):
    adj1 = set((genome1[i], genome1[i+1]) for i in range(len(genome1)-1))
    adj2 = set((genome2[i], genome2[i+1]) for i in range(len(genome2)-1))
    
    breakpoints = adj1.symmetric_difference(adj2)
    
    return len(breakpoints)

def find_breakpoints(genome1, genome2):
    adj1 = set((genome1[i], genome1[i+1]) for i in range(len(genome1)-1))
    adj2 = set((genome2[i], genome2[i+1]) for i in range(len(genome2)-1))
    
    breakpoints = adj1.symmetric_difference(adj2)
    
    print(f"\n✅ Total number of Breakpoints: {len(breakpoints)}\n")
    print("Breakpoints (gene pairs that are broken):")
    for b in breakpoints:
        print(b)

find_breakpoints(genome1, genome2)

# Focus on the area around the mutation for better visualization
buffer = 20  # Add some genes before and after mutation
vis_start = max(0, start - buffer)
vis_end = min(len(genome1), end + buffer)

small_genome1 = genome1[vis_start:vis_end]
small_genome2 = genome2[vis_start:vis_end]

def verify_mutation(genome1, genome2, start_index, end_index):
    print("\n=== MUTATION VERIFICATION ===")
    print(f"Original segment in genome1: {genome1[start_index:end_index]}")
    print(f"Inverted segment in genome2: {genome2[start_index:end_index]}")
    
    original = genome1[start_index:end_index]
    should_be = [-x for x in reversed(original)]
    actual = genome2[start_index:end_index]
    
    if should_be == actual:
        print("✅ Mutation correctly applied!")
    else:
        print("❌ Mutation NOT correctly applied!")
        print(f"Expected: {should_be}")
        print(f"Actual:   {actual}")
    
    before_mutation = genome1[start_index-1] if start_index > 0 else None
    after_mutation = genome1[end_index] if end_index < len(genome1) else None
    
    print("\nBreakpoints created by mutation:")
    if start_index > 0:
        print(f"Before inversion: {genome1[start_index-1]} -> {genome1[start_index]}")
        print(f"After inversion:  {genome2[start_index-1]} -> {genome2[start_index]}")
    
    if end_index < len(genome1):
        print(f"Before inversion: {genome1[end_index-1]} -> {genome1[end_index]}")
        print(f"After inversion:  {genome2[end_index-1]} -> {genome2[end_index]}")

# Verify the mutation first
verify_mutation(genome1, genome2, start, end)

def draw_breakpoint_graph(small_genome1, small_genome2):
    G = nx.Graph()
    
    # Create directed adjacencies with tuple order preserved (important!)
    adj1 = [(small_genome1[i], small_genome1[i+1]) for i in range(len(small_genome1)-1)]
    adj2 = [(small_genome2[i], small_genome2[i+1]) for i in range(len(small_genome2)-1)]
    
    # Store as strings to preserve direction (crucial fix)
    adj1_str = set(f"{u}->{v}" for u, v in adj1)
    adj2_str = set(f"{u}->{v}" for u, v in adj2)
    
    # Calculate unique edges
    only_in_genome1 = adj1_str - adj2_str
    only_in_genome2 = adj2_str - adj1_str
    in_both = adj1_str & adj2_str
    
    # Create a mapping from string representation back to tuples
    edges1 = {f"{u}->{v}": (u, v) for u, v in adj1}
    edges2 = {f"{u}->{v}": (u, v) for u, v in adj2}
    
    # Add all nodes
    all_nodes = set()
    for u, v in adj1 + adj2:
        all_nodes.add(u)
        all_nodes.add(v)
    
    for node in all_nodes:
        G.add_node(node)
    
    # Add edges with appropriate colors
    for edge_str in only_in_genome1:
        u, v = edges1[edge_str]
        G.add_edge(u, v, color='red')
    
    for edge_str in only_in_genome2:
        u, v = edges2[edge_str]
        G.add_edge(u, v, color='blue')
    
    for edge_str in in_both:
        u, v = edges1[edge_str] if edge_str in edges1 else edges2[edge_str]
        G.add_edge(u, v, color='purple')
    
    # Get edge colors for drawing
    colors = [G[u][v]['color'] for u, v in G.edges()]
    
    # Print debugging info
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

# Draw the graph AFTER verifying
draw_breakpoint_graph(small_genome1, small_genome2)
