import pandas as pd
import os
import random
import networkx as nx
import matplotlib.pyplot as plt

def load_genome_from_file(filepath):
    """
    Load a genome from a text file containing space-separated gene numbers
    """
    try:
        with open(filepath, 'r') as f:
            genome_text = f.read().strip()
            genome = [int(x) for x in genome_text.split()]
        return genome
    except Exception as e:
        print(f"Error loading genome from {filepath}: {e}")
        return None

def parse_gtf_file(gtf_file):
    """
    Parse a GTF file and return a list of signed genome elements
    """
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
    
    genome = []
    for idx, row in df.iterrows():
        gene_num = gene_to_num[row['gene_id']]
        if row['strand'] == '+':
            genome.append(gene_num)
        else:
            genome.append(-gene_num)
    
    return genome

def count_breakpoints(genome1, genome2):
    """
    Count the number of breakpoints between two genomes
    """
    adj1 = set((genome1[i], genome1[i+1]) for i in range(len(genome1)-1))
    adj2 = set((genome2[i], genome2[i+1]) for i in range(len(genome2)-1))
    
    breakpoints = adj1.symmetric_difference(adj2)
    
    return len(breakpoints)

def find_breakpoints(genome1, genome2):
    """
    Find and list the breakpoints between two genomes
    """
    adj1 = set((genome1[i], genome1[i+1]) for i in range(len(genome1)-1))
    adj2 = set((genome2[i], genome2[i+1]) for i in range(len(genome2)-1))
    
    breakpoints = adj1.symmetric_difference(adj2)
    
    print(f"\n✅ Total number of Breakpoints: {len(breakpoints)}\n")
    print("Breakpoints (gene pairs that are broken):")
    for b in breakpoints:
        print(b)
    
    return breakpoints

def simulate_inversion(genome1):
    """
    Create a new genome by simulating an inversion mutation
    """
    genome2 = genome1.copy()
    
    # Create a mutation
    start = random.randint(0, len(genome2) - 50)  # Start point
    end = start + random.randint(20, 30)          # Length of 20 to 30 genes
    
    # Invert the region
    genome2[start:end] = [-x for x in reversed(genome2[start:end])]
    
    print(f"\n✅ Simulated Mutation: Inverted genes from index {start} to {end-1}.")
    
    return genome2, start, end

def verify_mutation(genome1, genome2, start_index, end_index):
    """
    Verify that a mutation was correctly applied
    """
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
    
    print("\nBreakpoints created by mutation:")
    if start_index > 0:
        print(f"Before inversion: {genome1[start_index-1]} -> {genome1[start_index]}")
        print(f"After inversion:  {genome2[start_index-1]} -> {genome2[start_index]}")
    
    if end_index < len(genome1):
        print(f"Before inversion: {genome1[end_index-1]} -> {genome1[end_index]}")
        print(f"After inversion:  {genome2[end_index-1]} -> {genome2[end_index]}")

def draw_breakpoint_graph(genome1, genome2, focus_region=None):
    """
    Draw the breakpoint graph between two genomes
    
    Parameters:
    - genome1, genome2: The two genomes to compare
    - focus_region: Optional tuple (start, end) to focus on a specific region
    """
    # If focus_region is provided, use it to limit the genomes
    if focus_region:
        start, end = focus_region
        small_genome1 = genome1[start:end]
        small_genome2 = genome2[start:end]
    else:
        # Default: use first 200 genes or center on mutation if detected
        if len(genome1) > 400:
            # Try to find mutation region by checking for differences
            diff_indices = []
            for i in range(min(len(genome1), len(genome2))):
                if genome1[i] != genome2[i]:
                    diff_indices.append(i)
            
            if diff_indices:
                # Center on differences
                center = sum(diff_indices) // len(diff_indices)
                start = max(0, center - 50)
                end = min(len(genome1), center + 50)
                small_genome1 = genome1[start:end]
                small_genome2 = genome2[start:end]
            else:
                # No differences found, use first 200
                small_genome1 = genome1[:200]
                small_genome2 = genome2[:200]
        else:
            # Genomes are small enough to show whole
            small_genome1 = genome1
            small_genome2 = genome2
    
    G = nx.Graph()
    
    # Create directed adjacencies
    adj1 = [(small_genome1[i], small_genome1[i+1]) for i in range(len(small_genome1)-1)]
    adj2 = [(small_genome2[i], small_genome2[i+1]) for i in range(len(small_genome2)-1)]
    
    # Store as strings to preserve direction
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

def main():
    # Get current working folder
    current_folder = os.getcwd()
    
    # Process command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Compare genomes and generate breakpoint graphs')
    parser.add_argument('--gtf', type=str, default="gencode.v47.basic.annotation.gtf", 
                        help='GTF file to parse (default: gencode.v47.basic.annotation.gtf)')
    parser.add_argument('--genome2', type=str, default=None, 
                        help='Path to file containing second genome (optional)')
    args = parser.parse_args()
    
    gtf_file = os.path.join(current_folder, args.gtf)
    
    print(f"Parsing GTF file: {gtf_file}")
    genome1 = parse_gtf_file(gtf_file)
    
    print("\n✅ Signed Genome1 (first 100 genes):")
    print(genome1[:100])
    
    # Save Genome1 to file
    with open('signed_genome1_list.txt', 'w') as f:
        f.write(' '.join(map(str, genome1)))
    
    print("\nSigned Genome1 saved as 'signed_genome1_list.txt'.")
    
    # Either load genome2 from file or simulate it
    if args.genome2 and os.path.exists(args.genome2):
        print(f"\nLoading Genome2 from file: {args.genome2}")
        genome2 = load_genome_from_file(args.genome2)
        
        if genome2:
            print("\n✅ Signed Genome2 (first 100 genes):")
            print(genome2[:100])
            
            # Compare genomes
            find_breakpoints(genome1, genome2)
            
            # Draw breakpoint graph
            draw_breakpoint_graph(genome1, genome2)
        else:
            print("Failed to load genome2. Exiting.")
            return
    else:
        # Simulate a mutation
        genome2, start, end = simulate_inversion(genome1)
        
        print("\n✅ Signed Genome2 (first 100 genes):")
        print(genome2[:100])
        
        # Save Genome2 to file
        with open('signed_genome2_list.txt', 'w') as f:
            f.write(' '.join(map(str, genome2)))
        
        print("\nSigned Genome2 saved as 'signed_genome2_list.txt'.")
        
        # Verify the mutation
        verify_mutation(genome1, genome2, start, end)
        
        # Find breakpoints
        find_breakpoints(genome1, genome2)
        
        # Draw breakpoint graph with focus on mutation region
        buffer = 20
        vis_start = max(0, start - buffer)
        vis_end = min(len(genome1), end + buffer)
        draw_breakpoint_graph(genome1, genome2, (vis_start, vis_end))

if __name__ == "__main__":
    main()
