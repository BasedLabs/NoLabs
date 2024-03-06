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
      <div class="text-h7 q-pl-sm q-pb-md text-info" ><a class="text-light-blue" :href="this.ligand.link" target="_blank">{{ this.ligand.link }}</a></div>
      <div class="text-subtitle1">SMILES: {{ ligand.data.smiles }}</div>
      <PdbViewer v-if="sdfFile" :sdf-file="sdfFile" />
    </q-card-section>
  </q-card>
</template>

<script lang="ts">
import { QCard, QCardSection } from "quasar";

import PdbViewer from 'src/components/PdbViewer.vue'
import {defineComponent, PropType} from "vue";
import {ExtendedLigandMetaData} from "../targets/types";

export default defineComponent({
  name: "LigandDetail",
  components: {
    PdbViewer,
    QCard,
    QCardSection,
  },
  props: {
    originalLigand: {
      type: Object as PropType<ExtendedLigandMetaData>,
      required: true
    },
  },
  data() {
    return {
      ligand: this.originalLigand,
      sdfFile: null as File | null
    };
  },
  methods: {
  },
  mounted() {
    if (this.ligand.data?.sdf_file) {
      this.sdfFile= new File([new Blob([this.ligand.data.sdf_file])], "ligand.sdf");
    }
  },
})
</script>
