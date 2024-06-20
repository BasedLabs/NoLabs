<template>
    <q-card flat bordered v-if="!ligand">
      <q-card-section>
        <div class="text-caption">Smiles String</div>
      </q-card-section>
      <q-separator />

      <q-card-section>
        <q-skeleton />
      </q-card-section>
      <q-card-section>
        <div class="text-caption">3D View</div>
      </q-card-section>
      <q-skeleton height="200px" square />
      <q-separator />
    </q-card>
    <q-card v-if="ligand">
      <q-card-section>
        <div class="text-h6">Ligand: {{ ligand.name }}</div>
        <div class="text-h7 q-pl-sm q-pt-md q-pb-md text-info" ><a class="text-light-blue" :href="ligand.link" target="_blank">{{ ligand.link }}</a></div>
        <div class="text-subtitle1">SMILES: {{ ligand.smiles_content }}</div>
        <PdbViewer v-if="sdfFile" :sdf-file="sdfFile" :key="sdfFile.size" />
      </q-card-section>
    </q-card>
  </template>

  <script lang="ts">
  import { QCard, QCardSection } from "quasar";

  import PdbViewer from 'src/components/PdbViewer.vue'
  import {defineComponent, PropType} from "vue";
  import { LigandContentResponse } from "src/refinedApi/client";

  export default defineComponent({
    name: "LigandDetail",
    components: {
      PdbViewer,
      QCard,
      QCardSection,
    },
    props: {
      originalLigand: {
        type: Object as PropType<LigandContentResponse>,
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
      if (this.ligand.sdf_content) {
        this.sdfFile = new File([new Blob([this.ligand.sdf_content])], this.ligand.name + ".sdf");
      }
    },
  })
  </script>
