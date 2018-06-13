#! /usr/bin/perl

use strict;

use Getopt::Long;
use Bio::AlignIO;


my $compara_url;
my $super_trees = 0;
my $member_type = 'protein';
my $db_version = 85;
#my $disconnect_when_inactive = 0;

# Use -url mysql://anonymous@mysql.ebi.ac.uk:4157/ensembl_compara_fungi_15_68 to access EnsemblGenomes
GetOptions(
    'url=s'         => \$compara_url,
    'member_type=s' => \$member_type,
    'super_trees'   => \$super_trees,
    'db_version=i'   => \$db_version,
#    'disconnect_when_inactive'  => \$disconnect_when_inactive,
);

my $gene_tree_adaptor;

if ($compara_url) {

    use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
    my $dba = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(-url => $compara_url, -db_version=>$db_version);
    $gene_tree_adaptor = $dba->get_GeneTreeAdaptor;

} else {

    use Bio::EnsEMBL::Registry;
    my $reg = "Bio::EnsEMBL::Registry";
    $reg->load_registry_from_db(
            -host=>"ensembldb.ensembl.org",
            -user => "anonymous",
            -port => 5306,
            -db_version => $db_version,
            );
    $gene_tree_adaptor = $reg->get_adaptor("Multi", "compara", "GeneTree");

}
#$gene_tree_adaptor->disconnect_when_inactive($disconnect_when_inactive);


sub print_tree {
    my $node = shift;
    my $indent = shift;
    
    if ($node->is_leaf and not $node->isa('Bio::EnsEMBL::Compara::GeneTreeMember')) {
        # Special case for super-trees
        $gene_tree_adaptor->fetch_all_children_for_node($node);
    }
    if ($node->get_child_count == 1) {
        my $subtree = $gene_tree_adaptor->fetch_subtree_under_node($node->children->[0]);
        return print_tree($subtree, $indent);
    }

    print_node($node, $indent);

    $indent .= "\t";
    foreach my $child_node (@{$node->sorted_children}) {
        print $indent, "len\t", $child_node->distance_to_parent(), "\n";
        print_tree($child_node, $indent);
    }
}

my %cache_tn;
sub cached_taxon_name {
    my $node = shift;
    return $cache_tn{$node->taxon_id} if exists $cache_tn{$node->taxon_id};
    $cache_tn{$node->taxon_id} = $node->taxon->name;
    return $cache_tn{$node->taxon_id};
}

sub print_node {

    my $node = shift;
    my $indent = shift;

    print $indent."id\t".($node->node_id)."\n";
    print $indent."info\t{";

    my $kindofDup = $node->get_tagvalue("node_type");
    my $dup = 0;
    $dup = 1 if $kindofDup eq 'dubious';
    $dup = 2 if $kindofDup eq 'duplication';
    $dup = 10 if $kindofDup eq 'gene_split';
    print sprintf("'Duplication': %d, ", $dup);
    if ($dup) {
        print sprintf("'duplication_confidence_score': %f, ", $node->get_tagvalue("duplication_confidence_score") || 0);
    }

    if ($node->isa('Bio::EnsEMBL::Compara::GeneTreeMember')) {
        # leaf
        print sprintf("'protein_name': '%s', ", $node->name);
        print sprintf("'gene_name': '%s', ", $node->gene_member->stable_id);
        print sprintf("'taxon_id': %d, 'taxon_name': '%s'", $node->taxon_id, cached_taxon_name($node));
    } else {
        # internal node
        print sprintf("'tree_name': '%s', ", $node->tree->stable_id) if ($node->tree->root eq $node) and (defined $node->tree->stable_id);
        print sprintf("'Bootstrap': %d, ", $node->get_tagvalue("bootstrap")) if $node->has_tag("bootstrap");
        print sprintf("'taxon_id': %d, 'taxon_name': '%s'", $node->get_tagvalue("taxon_id"), $node->get_tagvalue("taxon_name"));
    }

    print "}\n";
}


if ($super_trees) {
    my $allclustersets = $gene_tree_adaptor->fetch_all(-clusterset_id => 'default', -tree_type => 'clusterset', -member_type => $member_type);
    print STDERR scalar(@$allclustersets), " clustersets to dump\n";
    foreach my $clusterset (@$allclustersets) {
        my $alltrees = $gene_tree_adaptor->fetch_subtrees($clusterset);
        print STDERR scalar(@$alltrees), " trees to dump\n";
        foreach my $tree (@$alltrees) {
            if ($tree->tree_type eq 'supertree') {
                $tree->expand_subtrees();
                $tree->stable_id(sprintf('SUPERTREE%d', $tree->root_id));
            } else {
                $tree->preload();
            }
            print_tree($tree->root, '');

        $tree->release_tree();
        }
    }
} else {
    my $alltrees = $gene_tree_adaptor->fetch_all(-clusterset_id => 'default', -tree_type => 'tree', -member_type => $member_type);
    print STDERR scalar(@$alltrees), " trees to dump\n";
    foreach my $tree (@$alltrees) {
        $tree->preload(); # This seems to stay in memory after. How to clear it?
        print_tree($tree->root, '');
        $tree->release_tree();
    }
}

