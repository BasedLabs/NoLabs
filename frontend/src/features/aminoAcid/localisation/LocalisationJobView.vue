<template>
    <div v-if="job">
      <q-separator></q-separator>
      <q-layout container style="height: 100vh">
        <JobHeader :job-name="job!.name" :on-job-name-change-submit="onJobNameChange">
          <q-btn color="info" size="md" outline label="Localisation parameters"
                 @click="showInferenceForm = !showInferenceForm"/>
        </JobHeader>
        <q-page-container>
          <div class="row" v-if="jobHasGeneratedData">
            <div class="col-6">
              <div class="q-ma-sm">
                <AminoAcidTable :on-amino-acid-open="setActiveAminoAcid" :rows="job?.aminoAcids"/>
              </div>
            </div>
            <div class="col-2">
              <div class="q-ma-sm">
                <q-toolbar class="bg-primary text-white shadow-2">
                  <q-toolbar-title>Probabilities</q-toolbar-title>
                </q-toolbar>
                <q-list bordered separator>
                  <q-item clickable v-ripple v-for="data in listItemData"
                          @mouseenter="setActiveListItem(data.key)"
                          @mouseleave="resetActiveListItem" :key="data.key">
                    <q-item-section>
                      {{ data.text }}
                    </q-item-section>
                  </q-item>
                </q-list>
              </div>
            </div>
            <div class="col-4">
              <div class="q-pl-sm q-ma-sm">
                <q-img
                    :src="activeImage"
                    style="height: 500px; max-width: 500px;"
                />
              </div>
            </div>
          </div>
        </q-page-container>
      </q-layout>
      <q-dialog v-model="showInferenceForm" position="right" :maximized="true">
        <q-card>
          <q-card-section class="row items-center q-pb-none">
            <div class="text-h6">Localisation parameters</div>
            <q-space/>
            <q-btn icon="close" flat round dense v-close-popup/>
          </q-card-section>
          <q-card-section>
            <AminoAcidInferenceForm :on-submit="onSubmit" :properties="job?.properties"/>
          </q-card-section>
        </q-card>
      </q-dialog>
    </div>
  </template>

  <script lang="ts">
  import {defineComponent, ref} from 'vue'
  import {QVueGlobals, QSpinnerOrbit} from 'quasar';
  import useLocalisationStore from "src/features/aminoAcid/localisation/storage";
  import {AminoAcid} from "src/features/aminoAcid/localisation/types";
  import {Job} from "src/features/aminoAcid/types";
  import AminoAcidInferenceForm from "src/features/aminoAcid/AminoAcidInferenceForm.vue";
  import AminoAcidTable from "src/features/aminoAcid/AminoAcidTable.vue";
  import JobHeader from "../../../components/JobHeader.vue";


  export default defineComponent({
    name: 'LocalisationJobView',
    data() {
      const store = useLocalisationStore();

      const images: { [key: string]: string } = {
        mitochondria: '/localisationImages/mithochondria.png',
        nucleus: '/localisationImages/nucleus.png',
        cytoplasm: '/localisationImages/cytoplasm.png',
        extracellular: '/localisationImages/original.png',
        other: '/localisationImages/original.png',
        original: '/localisationImages/original.png',
      };

      return {
        job: null as Job<AminoAcid>,
        showInferenceForm: false,
        store,
        activeImage: images.original,
        listItemData: [] as Array<{
          key: string,
          text: string,
          image: string
        }>,
        images
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
      setActiveAminoAcid(aminoAcidName: string): void {
        const aminoAcid = this.job?.aminoAcids.find(x => x.name === aminoAcidName);
        this.listItemData = [
          {
            key: 'mitochondria',
            text: `Mitochondria ${this.formatText(aminoAcid?.mitochondialProteins!)}`,
            image: this.images.mitochondria
          },
          {
            key: 'nucleus',
            text: `Nucleus  ${this.formatText(aminoAcid?.nuclearProteins!)}`,
            image: this.images.nucleus
          },
          {
            key: 'cytoplasm',
            text: `Cytoplasm ${this.formatText(aminoAcid?.cytosolicProteins!)}`,
            image: this.images.cytoplasm
          },
          {
            key: 'extracellular',
            text: `Extracellular ${this.formatText(aminoAcid?.extracellularSecretedProteins!)}`,
            image: this.images.original
          },
          {
            key: 'other',
            text: `Other ${this.formatText(aminoAcid?.otherProteins!)}`,
            image: this.images.other
          }
        ]
        this.activeImage = this.images.original;
      },
      setJob(job: Job<AminoAcid> | null){
        if (job !== null) {
          this.job = job;

          if(job.aminoAcids.length > 0){
            this.setActiveAminoAcid(job.aminoAcids[0]!.name);
          }
        }
      },
      setActiveListItem(key: string) {
        this.activeImage = this.listItemData.find(x => x.key === key)!.image;
      },
      resetActiveListItem() {
        this.activeImage = this.images.original;
      },
      formatText(prob: number): number {
        return (prob as any).toFixed(6) * 100;
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
          fastas: data.fastas
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
      JobHeader,
      AminoAcidTable,
      AminoAcidInferenceForm
    }
  })
  </script>
