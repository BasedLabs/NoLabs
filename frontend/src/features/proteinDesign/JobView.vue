<template>
  <div v-if="jobLoaded">
    <q-separator></q-separator>
    <q-layout container style="height: 100vh">
      <JobHeader :job-name="job!.name" :on-job-name-change-submit="onJobNameChange">
        <q-btn color="info" size="md" outline label="Binder parameters"
               @click="showInferenceForm = !showInferenceForm"/>
      </JobHeader>
      <q-page-container>
        <div class="row" v-if="jobHasGeneratedData">
          <div class="col-5">
            <div class="q-ma-sm">
              <PdbViewer :pdb-file="job?.properties.inputPdbFile"/>
            </div>
          </div>
          <div class="col-2">
            <div class="q-pl-sm q-ma-sm">
              <q-table
                  title="Generated pdbs"
                  :rows="generatedPdbsTableRows"
                  :columns="generatedPdbsTableColumns"
                  row-key="name"
                  @row-click="(_, row) => {generatedPdbsReload = !generatedPdbsReload; selectedGeneratedPdbIndex = row.id;}"
                  :pagination="{rowsPerPage: 5}"
              />
            </div>
          </div>
          <div class="col-5">
            <div class="q-mt-sm q-mb-sm q-mr-sm">
              <PdbViewer :pdb-file="job!.generatedPdbs[selectedGeneratedPdbIndex]"
                         :key="generatedPdbsReload"/>
            </div>
          </div>
        </div>
      </q-page-container>
    </q-layout>
    <q-dialog v-model="showInferenceForm" position="right" :maximized="true">
      <q-card>
        <q-card-section class="row items-center q-pb-none">
          <div class="text-h6">Protein binder design parameters</div>
          <q-space/>
          <q-btn icon="close" flat round dense v-close-popup/>
        </q-card-section>
        <q-card-section>
          <InferenceFormView :on-submit="onSubmit" :properties="job!.properties"/>
        </q-card-section>
      </q-card>
    </q-dialog>
  </div>
</template>

<script lang="ts">
import {defineComponent, ref} from 'vue'
import useProteinDesignStore from 'src/features/proteinDesign/storage';
import {QVueGlobals, QSpinnerOrbit} from 'quasar';
import InferenceFormView from "src/features/proteinDesign/InferenceFormView.vue";
import {Job, JobProperties} from "src/features/proteinDesign/types";
import PdbViewer from "src/components/PdbViewer.vue";
import JobHeader from "src/components/JobHeader.vue";


export default defineComponent({
  name: 'ProteinDesignJobView',
  props: {
    experimentId: {
      type: String,
      required: true,
    },
  },
  data() {
    const store = useProteinDesignStore();

    return {
      job: null as Job,
      showInferenceForm: false,
      store,
      selectedGeneratedPdbIndex: 0,
      generatedPdbsReload: false
    }
  },
  computed: {
    generatedPdbsTableColumns() {
      return [
        {
          name: 'id',
          required: true,
          label: '#',
          align: 'left',
          sortable: false,
          field: row => row.id,
        },
        {
          name: 'name',
          required: true,
          label: 'Name',
          align: 'left',
          sortable: false,
          field: row => row.name
        }
      ]
    },
    generatedPdbsTableRows(): { id: number, name: string }[] {
      return this.job ? this.job.generatedPdbs.map((file, i) => {
        return {
          id: i, name: file.name
        }
      }) : [];
    },
    jobLoaded(): boolean {
      return this.job !== null;
    },
    jobHasGeneratedData(): boolean {
      return this.jobLoaded && this.job!.generatedPdbs.length > 0;
    }
  },
  methods: {
    async onJobNameChange(newJobName: string) {
      await this.store.changeJobName(this.job?.id as string, newJobName);
      this.job!.name = newJobName;
    },
    async onSubmit(properties: JobProperties) {
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: 'Running AI models. This can take a couple of minutes'
      });

      const response = await this.store.inference({
        experimentId: this.experimentId,
        jobId: this.job?.id,
        jobName: this.job?.name as string,
        pdbFile: properties.inputPdbFile!,
        contig: properties.contig,
        numberOfDesigns: properties.numberOfDesigns,
        timesteps: properties.timesteps,
        hotspots: properties.hotspots
      });

      if (response.job !== null) {
        this.job = response.job;
      }

      this.showInferenceForm = false;

      this.$q.loading.hide();
    },
    changeJobName() {
      this.$q.dialog({
        color: 'info',
        title: 'Prompt',
        message: 'Enter new job name',
        prompt: {
          model: this.job!.name,
          required: true,
          type: 'text' // optional
        },
        cancel: true,
        persistent: true
      }).onOk(async data => {
        if (!data)
          return;
        this.$q.loading.show({
          spinner: QSpinnerOrbit,
          message: 'Changing job name'
        });
        await this.store.changeJobName(this.job?.id as string, data);
        this.job!.name = data;
        this.$q.loading.hide();
      });
    }
  },
  async mounted() {
    const jobId = this.$route.params.jobId as string;

    this.$q.loading.show({
      spinner: QSpinnerOrbit,
      message: `Job ${jobId}`
    });

    const response = await this.store.getJob(jobId);

    if (response.job !== null) {
      this.job = response.job;
    }

    this.$q.loading.hide();

    if (!this.jobHasGeneratedData) {
      this.showInferenceForm = true;
    }
  },
  components: {
    JobHeader,
    PdbViewer,
    InferenceFormView
  }
})
</script>
