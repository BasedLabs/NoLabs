<script lang="ts">
import {defineComponent, PropType} from "vue";
import PdbViewer from '../../../components/PdbViewer.vue';
import useSmallMoleculesDesignStore from "../storage";

interface Data {
  smilesData: Smiles[];
  params: Params | null;
  interval: any;
}

export default defineComponent({
  name: "TableJobExpandedRow",
  components: {PdbViewer},
  store: useSmallMoleculesDesignStore(),
  smilesColumns: [
    {
      name: 'smiles',
      label: 'Formula',
      align: 'left',
      sortable: false
    },
    {
      name: 'drugLikeness',
      label: 'Drug Likeness',
      align: 'left',
      sortable: true
    },
    {
      name: 'score',
      label: 'Score',
      align: 'left',
      sortable: true
    },
    {
      name: 'createdAt',
      label: 'Generated At',
      align: 'left',
      sortable: true
    }
  ],
  props: {
    job: {
      type: Object as PropType<Job>,
      required: true
    },
    expanded: {
      type: Boolean,
      required: true
    },
    onSubmit: {
      type: Function as PropType<(params: Params) => Promise<void>>,
      required: true
    }
  },
  data(): Data {
    return {
      smilesData: [] as Smiles[],
      params: null as Params | null,
      interval: null
    }
  },
  computed: {
    readonly() {
      return !this.job.running && !this.job.learningCompleted;
    }
  },
  methods: {
    async submit() {
      await this.onSubmit(this.params!);
    }
  },
  mounted() {
    this.interval = setInterval(async () => {
      if(this.expanded){
        this.smilesData = this.$options.store.smilesData(this.job.id);
      }
    }, 1000);
  },
  unmounted() {
    if(this.interval){
      clearInterval(this.interval);
    }
  }
})
</script>

<template>
  <div class="row">
    <div class="col-sm-12 col-md-4">
      <div class="q-pa-md">
        <q-table
          title="Generated smiles"
          :rows="smilesData"
          :columns="$options.smilesColumns"
          row-key="name"
        />
      </div>
    </div>
    <div class="col-sm-12 col-md-5">
      <PdbViewer v-if="params?.pdbFile" :pdb-file="params?.pdbFile"/>
    </div>
    <div class="col-sm-12 col-md-3">
      <q-form v-if="params" class="q-gutter-md" @submit="submit">
        <q-file v-if="params" filled bottom-slots accept=".fasta" v-model="params.pdbFile"
                :readonly="readonly"
                label=".pdb file" counter>
          <template v-slot:prepend>
            <q-icon name="cloud_upload" @click.stop.prevent/>
          </template>
          <template v-slot:append>
            <q-icon name="close" class="cursor-pointer"/>
          </template>
          <template v-slot:hint>
            Upload .pdb file
          </template>
        </q-file>
        <p>Search space</p>
        <q-input filled v-model="params.centerX" label="Center X" type="number" step="1" :readonly="readonly"
                 :rules="[val => val && val > 0 || 'Please type something']">
        </q-input>
        <q-input filled v-model="params.centerY" label="Center Y" type="number" step="1" :readonly="readonly"
                 :rules="[val => val && val > 0 || 'Please type something']">
        </q-input>
        <q-input filled v-model="params.centerZ" label="Center Z" type="number" step="1" :readonly="readonly"
                 :rules="[val => val && val > 0 || 'Please type something']">
        </q-input>
        <q-input filled v-model="params.centerX" label="Size X" type="number" step="1" :readonly="readonly"
                 :rules="[val => val && val > 0 || 'Please type something']">
        </q-input>
        <q-input filled v-model="params.centerY" label="Size Y" type="number" step="1" :readonly="readonly"
                 :rules="[val => val && val > 0 || 'Please type something']">
        </q-input>
        <q-input filled v-model="params.centerZ" label="Size Z" type="number" step="1" :readonly="readonly"
                 :rules="[val => val && val > 0 || 'Please type something']">
        </q-input>
        <p>Other parameters</p>
        <p>1 epoch training will take 20-30 min approx</p>
        <q-input filled v-model="params.batchSize" label="Epoch batch size" type="number" step="1" :readonly="readonly"
                 :rules="[val => val && val > 0 || 'Please type something']">
        </q-input>
        <q-input filled v-model="params.minscore" label="Minimum docking acceptance score" type="number" step="1" :readonly="readonly"
                 :rules="[val => val && val > 0 || 'Please type something']">
        </q-input>
        <q-input filled v-model="params.epochs" label="Epochs" type="number" step="1" :readonly="readonly"
                 :rules="[val => val && val > 0 || 'Please type something']">
        </q-input>
        <div>
          <q-btn label="Save" size="large" type="submit" color="info" :disabled="readonly"/>
        </div>
      </q-form>
    </div>
  </div>
</template>
