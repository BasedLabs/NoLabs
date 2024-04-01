<template>
  <div class="q-pa-md">
    <q-table
      flat bordered
      title="Design jobs"
      :rows="jobs"
      :columns="columns"
      row-key="name"
      separator="none"
    >

      <template v-slot:header="props">
        <q-tr :props="props">
          <q-th auto-width/>
          <q-th
            v-for="col in props.cols"
            :key="col.name"
            :props="props"
          >
            {{ col.label }}
          </q-th>
        </q-tr>
      </template>

      <template v-slot:body="props">
        <q-tr :props="props">
          <q-td auto-width>
            <q-btn size="sm" color="info" square outline dense @click="props.expand = !props.expand"
                   :icon="props.expand ? 'remove' : 'add'"/>
          </q-td>
          <q-td
            v-for="col in props.cols"
            :key="col.name"
            :props="props"
          >
            {{ col.value }}
          </q-td>
          <q-td auto-width>
            <q-btn size="md" color="info" square outline dense @click="openLogs(props.row.id)">Logs
            </q-btn>
          </q-td>
          <JobControlButtons :job="props.row" :resume-job="generateMoreMolecules" :stop-job="stopJob"
                             :start-job="startJob"/>
          <q-td auto-width class="q-ml-md">
            <q-btn size="md" color="negative" square dense icon="delete"/>
          </q-td>
          <q-td auto-width class="q-ml-md" v-if="props.row.running">
            Running
            <q-spinner-oval
              color="info"
              size="3em"
            />
          </q-td>
        </q-tr>
        <q-tr v-show="props.expand" :props="props">
          <q-td colspan="100%">
            <div class="text-left">This is expand slot for row above: {{ props.row.name }}.</div>
          </q-td>
        </q-tr>
      </template>
    </q-table>
    <LogsModal :job-id="logsModalJobId" v-model:visible="logsModalVisible"/>
  </div>
</template>


<script lang="ts">
import {defineComponent} from "vue";
import useSmallMoleculesDesignStore from "./storage";
import {QSpinnerOrbit} from "quasar";
import LogsModal from "./components/LogsModal.vue";
import JobControlButtons from "./components/JobControlButtons.vue";

const columns = [
  {name: 'id', align: 'center', label: 'ID', field: 'id', sortable: false},
  {name: 'name', align: 'center', label: 'Name', field: 'name', sortable: false},
  {name: 'createdAt', align: 'center', label: 'Created At', field: 'createdAt', sortable: true}
]

interface Data {
  jobs: Job[]
  columns: any;
  logsModalVisible: boolean;
  logsModalJobId: string | null;
  logsContent: Logs | null;
  interval: any;
}

export default defineComponent({
  components: {JobControlButtons, LogsModal},
  store: useSmallMoleculesDesignStore(),
  name: 'SmallMoleculesDesignExperiment',
  data(): Data {
    return {
      columns,
      jobs: [] as Job[],
      logsModalVisible: false,
      logsModalJobId: null,
      logsContent: null as Logs | null,
      interval: null
    }
  },
  methods: {
    openLogs(id: string) {
      this.logsModalJobId = id;
      this.logsModalVisible = !this.logsModalVisible;
    },
    async loader(title: string, operation: any) {
      try {
        this.$q.loading.show({
          spinner: QSpinnerOrbit,
          message: title
        });
        await operation();
      } finally {
        this.$q.loading.hide();
      }
    },
    async startJob(id: string) {
      await this.loader('Starting job', async () => {
        await this.$options.store.startJob(id);
        const job = this.jobs.find(x => x.id === id);
        const index = this.jobs.indexOf(job!);
        this.jobs[index] = await this.$options.store.job(id);
      });
    },
    async stopJob(id: string) {
      this.$q.dialog({
        title: 'Confirm',
        message: 'Would you like to stop the job? The progress will be lost.',
        cancel: true,
        persistent: true
      }).onOk(async () => {
        await this.loader('Stopping job', async () => {
          await this.$options.store.stopJob(id);
          const job = this.jobs.find(x => x.id === id);
          const index = this.jobs.indexOf(job!);
          this.jobs[index] = await this.$options.store.job(id);
        });
      });
    },
    async generateMoreMolecules(id: string) {
      await this.loader('Starting job', async () => {
        await this.$options.store.generateMoreMolecules(id);
        const job = this.jobs.find(x => x.id === id);
        const index = this.jobs.indexOf(job!);
        this.jobs[index] = await this.$options.store.job(id);
      });
    },
    deleteJob(id: string) {
      this.$q.dialog({
        title: 'Delete?',
        cancel: true,
        persistent: true
      }).onOk(async () => {
        await this.loader('Deleting job', async () => {
          await this.$options.store.deleteJob(id);
          const job = this.jobs.find(x => x.id === id);
          const index = this.jobs.indexOf(job!);
          this.jobs.splice(index, 1);
        });
      });
    },
    async loadJobs() {
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: 'Loading jobs'
      });

      this.jobs = await this.$options.store.jobs();

      this.$q.loading.hide();
    }
  },
  mounted() {
    this.interval = setInterval(async () => {
      await this.loadJobs();
    }, 1000);
  },
  unmounted() {
    clearInterval(this.interval);
  }
})
</script>
