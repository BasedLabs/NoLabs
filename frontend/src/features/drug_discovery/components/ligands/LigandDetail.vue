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
      <div ref="viewerContainer" class="ligand-viewer"></div>
    </q-card-section>
  </q-card>
</template>

<script>
import { QCard, QCardSection } from "quasar";
import * as NGL from "ngl";

export default {
  name: "LigandDetail",
  components: {
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
    createViewer() {
      if (this.viewer || !this.ligand.data.sdf_file) {
        return;
      }
      setTimeout(async () => {
        this.viewer = new NGL.Stage(this.$refs.viewerContainer);
        const stringBlob = new Blob([this.ligand.data.sdf_file], {
          type: "text/plain",
        });
        this.viewer.loadFile(stringBlob, { ext: "sdf" }).then((component) => {
          component.addRepresentation("ball+stick");
          component.autoView();
        });
      }, 100);
    },
  },
  mounted() {
    this.createViewer();
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
