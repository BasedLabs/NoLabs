<template>
  <div v-if="jobLoaded">
    <q-separator></q-separator>
    <q-layout container style="height: 100vh">
      <JobHeader :job-name="job!.name" :on-job-name-change-submit="onJobNameChange">
        <q-btn color="info" size="md" outline label="Gene ontology parameters"
               @click="showInferenceForm = !showInferenceForm"/>
      </JobHeader>
      <q-page-container>
        <div class="row" v-if="jobHasGeneratedData">
          <div :class="tiles.one.current"
               style="transition: all .1s linear;">
            <div class="q-ma-sm">
              <q-btn size="xs" flat color="info" style="width: 100%" class="q-mb-xs" label="EXPAND"
                     @click="expandTile('one');"/>
              <AminoAcidTable :on-amino-acid-open="setActiveAminoAcid" :rows="job!.aminoAcids"/>
            </div>
          </div>
          <div :class="tiles.two.current"
               style="transition: all .1s linear;">
            <div class="q-ma-sm">
              <q-btn size="xs" flat color="info" style="width: 100%" class="q-mb-xs" label="EXPAND"
                     @click="expandTile('two');"/>
              <GeneOntologyTree :obo-graph="activeAminoAcid!.go" :key="activeAminoAcid?.name"/>
            </div>
          </div>
        </div>
      </q-page-container>
    </q-layout>
    <q-dialog v-model="showInferenceForm" position="right" :maximized="true">
      <q-card>
        <q-card-section class="row items-center q-pb-none">
          <div class="text-h6">Gene ontology parameters</div>
          <q-space/>
          <q-btn icon="close" flat round dense v-close-popup/>
        </q-card-section>
        <q-card-section>
          <AminoAcidInferenceForm :on-submit="onSubmit" :properties="job!.properties"/>
        </q-card-section>
      </q-card>
    </q-dialog>
  </div>
</template>

<script lang="ts">
import {defineComponent, PropType} from 'vue'
import {QSpinnerOrbit, QVueGlobals} from 'quasar';
import JobHeader from "src/components/JobHeader.vue";
import {Job} from "src/features/aminoAcid/types";
import AminoAcidInferenceForm from "src/features/aminoAcid/AminoAcidInferenceForm.vue";
import GeneOntologyTree from "src/features/aminoAcid/geneOntology/GeneOntologyTree.vue";
import {AminoAcid} from "src/features/aminoAcid/geneOntology/types";
import useGeneOntologyStore from "src/features/aminoAcid/geneOntology/storage";
import AminoAcidTable from "src/features/aminoAcid/AminoAcidTable.vue";


export default defineComponent({
  name: 'GeneOntologyJobView',
  data() {
    const store = useGeneOntologyStore();
    return {
      tiles: {
        one: {
          current: 'col-md-3',
          hover: 'col-md-6',
          leave: 'col-md-3',
          otherHover: 'col-md-1'
        },
        two: {
          current: 'col-md-9',
          hover: 'col-md-11',
          leave: 'col-md-9',
          otherHover: 'col-md-6'
        }
      } as { [index: string]: { current: String, hover: String, leave: String, otherHover: String } },
      job: null as Job<AminoAcid>,
      showInferenceForm: false,
      store,
      activeAminoAcid: null as AminoAcid | null | undefined
    }
  },
  computed: {
    jobLoaded(): boolean {
      return this.job !== null;
    },
    jobHasGeneratedData(): boolean {
      return this.jobLoaded && this.job!.aminoAcids.length > 0;
    },
    aminoAcidRows() {
      return this.job?.aminoAcids;
    },
  },
  methods: {
    expandTile(index: string) {
      if (this.tiles[index].current === this.tiles[index].hover) {
        for (const key in this.tiles) {
          this.tiles[key].current = this.tiles[key].leave;
        }
        return;
      }

      for (const key in this.tiles) {
        if (key == index) {
          this.tiles[key].current = this.tiles[key].hover;
        } else {
          this.tiles[key].current = this.tiles[key].otherHover;
        }
      }
    },
    setActiveAminoAcid(aminoAcidName: string): void {
      this.activeAminoAcid = this.job?.aminoAcids.find(x => x.name === aminoAcidName);
    },
    setJob(job: Job<AminoAcid> | null) {
      if (job !== null) {
        this.job = job;

        if (job.aminoAcids.length > 0) {
          this.setActiveAminoAcid(job.aminoAcids[0]!.name);
        }
      }
    },
    async onJobNameChange(newJobName: string) {
      await this.store.changeJobName(this.job?.id as string, newJobName);
      this.job!.name = newJobName;
    },
    async onSubmit(data: { fastas: Array<File> }) {
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: 'Running AI models. This can take a couple of minutes'
      });

      const response = await this.store.inference({
        jobId: this.job!.id,
        jobName: this.job!.name,
        fastas: data.fastas,
      });

      this.setJob(response.job);

      this.showInferenceForm = false;

      this.$q.loading.hide();
    },
  },
  async mounted() {
    const jobId = this.$route.params.jobId as string;

    this.$q.loading.show({
      spinner: QSpinnerOrbit,
      message: `Job ${jobId}`
    });

    const response = await this.store.getJob(jobId);

    this.setJob(response.job);

    this.$q.loading.hide();

    if (!this.jobHasGeneratedData) {
      this.showInferenceForm = true;
    }
  },
  components: {
    AminoAcidTable,
    GeneOntologyTree,
    JobHeader,
    AminoAcidInferenceForm
  }
})
</script>
