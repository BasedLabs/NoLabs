<template>
  <q-card flat bordered v-if="!ligand">
    <q-card-section>
      <div class="text-caption">Smiles String</div>
    </q-card-section>
    <q-separator />

    <q-card-section>
      <q-skeleton :type="type" />
    </q-card-section>
    <q-card-section>
      <div class="text-caption">3D View</div>
    </q-card-section>
    <q-skeleton height="200px" square />
    <q-separator />
  </q-card>
  <q-card v-if="ligand">
    <q-card-section>
      <div class="text-h6">Ligand: {{ ligand.ligand_name }}</div>
      <div class="text-subtitle1">SMILES: {{ ligand.data.smiles }}</div>
      <PdbViewer v-if="ligand.data.sdf_file" :sdf-file="ligand.data.sdf_file" />
    </q-card-section>
  </q-card>
</template>

<script>
import { QCard, QCardSection } from "quasar";
import * as NGL from "ngl";

import PdbViewer from 'src/components/PdbViewer.vue'

export default {
  name: "LigandDetail",
  components: {
    PdbViewer,
    QCard,
    QCardSection,
  },
  props: {
    originalLigand: {
      type: Object,
      required: true
    },
  },
  data() {
    return {
      viewer: null,
      ligand: this.originalLigand
    };
  },
  methods: {
  },
  mounted() {
    this.ligand.data.sdf_file = new File([new Blob([this.ligand.data.sdf_file])], "ligand.sdf");
  },
  beforeUnmount() {
    if (this.viewer) {
      this.viewer.dispose();
    }
  },
};
</script>

<style>
.ligand-viewer {
  width: 100%;
  height: 400px;
}
</style>
