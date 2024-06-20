<template>
  <q-card flat bordered v-if="!ligand.data">
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
  <q-card v-if="ligand.data">
    <q-card-section>
      <div class="text-h6">Ligand: {{ ligand.metaData.name }}</div>
      <div class="text-h7 q-pl-sm q-pt-md q-pb-md text-info"><a class="text-light-blue" :href="ligand.data.link" target="_blank">{{ ligand.data.link }}</a></div>
      <div class="text-subtitle1">SMILES: {{ ligand.data.smiles_content }}</div>
      <PdbViewer v-if="sdfFile" :sdf-file="sdfFile" :key="sdfFile.size" />
    </q-card-section>
  </q-card>
</template>

<script lang="ts">
import { QCard, QCardSection } from "quasar";
import PdbViewer from "src/components/PdbViewer.vue";
import { LigandContentResponse } from "src/refinedApi/client";
import { defineComponent, PropType } from "vue";
import {getLigandContent} from "src/features/drug_discovery/refinedApi"
import { useWorkflowStore } from "src/features/drug_discovery/components/workflow/storage";

export default defineComponent({
  name: "LigandDetail",
  components: {
    PdbViewer,
    QCard,
    QCardSection,
  },
  props: {
    ligandId: {
      type: String,
      required: true
    },
  },
  data() {
    return {
      ligand: {
        metaData: {
          id: this.ligandId,
          name: '',
        },
        data: {
          smiles_content: '',
          sdf_content: '',
          link: ''
        }
      },
      sdfFile: null as File | null
    };
  },
  async mounted() {
    await this.loadLigandData();
  },
  methods: {
    async loadLigandData() {
      try {
        const ligandContent: LigandContentResponse = await getLigandContent(this.ligandId) as LigandContentResponse;
        this.ligand.metaData.name = ligandContent.name;
        this.ligand.data.smiles_content = ligandContent.smiles_content || '';
        this.ligand.data.sdf_content = ligandContent.sdf_content || '';
        this.ligand.data.link = ligandContent.link || '';

        if (this.ligand.data.sdf_content) {
          this.sdfFile = new File([new Blob([this.ligand.data.sdf_content])], this.ligand.metaData.name + ".sdf");
        }
      } catch (error) {
        console.error('Error loading ligand data:', error);
      }
    }
  }
})
</script>
