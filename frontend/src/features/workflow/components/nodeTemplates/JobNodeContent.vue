<template>
  <q-card>
    <NodeContentHeader :name="name" :description="description" :inputs="['FASTA']" :outputs="['PDB']"
      :imageSrc="imageSrc" />
    <NodeContentBody :nodeId="nodeId" :name="name" />
  </q-card>
</template>

<script lang="ts">
import { defineComponent, computed } from 'vue';
import NodeContentHeader from './NodeContentHeader.vue';
import NodeContentBody from './NodeContentBody.vue';

export default defineComponent({
  name: 'JobNodeContent',
  components: {
    NodeContentHeader,
    NodeContentBody,
  },
  props: {
    nodeId: {
      type: String,
      required: true
    },
    name: {
      type: String,
      required: true
    },
    description: {
      type: String,
      required: true
    }
  },
  computed: {
    imageSrc(): string {
      switch (this.name) {
        case 'Esmfold light' || 'Esmfold' || 'Rosettafold':
          return '/folding_pic.png';
        case 'DiffDock':
          return '/docking_pic.png';
        case 'Small molecules design':
          return '/small_molecule.png';
        case 'Protein binder design':
          return '/rf_diffusion_pic.png';
        case 'Localisation':
          return '/localisation_pic.png';
        case 'GeneOntology':
          return '/gene_ontology_pic.png';
        case 'Conformations':
          return '/conformations_pic.png';
        case 'Binding pockets':
          return '/pocket_prediction_pic.png';
        case 'Msa generation':
          return '/msa_pic.png';
        default:
          return '/protein_folding.png';
      }
    }
  },
  mounted() {
    console.log(this.name)
  }
});
</script>
