<script lang="ts">
import {defineComponent, PropType} from 'vue'
import PdbViewer from '../../../components/PdbViewer.vue';
import useSmallMoleculesDesignStore from "../storage";
import {QSpinnerOrbit} from "quasar";

interface Data {
  params: Params | null,
}

export default defineComponent({
  name: "RowExpandedContent",
  components: {PdbViewer},
  store: useSmallMoleculesDesignStore(),
  props: {
    job: {
      type: Object as PropType<Job>,
      required: true
    }
  },
  data(): Data {
    return {
      params: null as Params | null,
    }
  },
  computed: {
    changeable(): boolean {
      return !!this.job && !this.job.running;
    }
  },
  async mounted() {
    this.$q.loading.show({
      spinner: QSpinnerOrbit,
      message: 'Loading data'
    });

    this.params = this.$options.store.params(this.job.id);

    this.$q.loading.hide();
  }
})
</script>

<template>
  <div class="row">
    <div class="col-sm-12 col-md-5"></div>
    <div class="col-sm-12 col-md-5">
      <q-file v-if="changeable" filled multiple bottom-slots accept=".pdb" v-model="params!.pdbFile"
              label=".pdb file" counter>
        <template v-slot:prepend>
          <q-icon name="cloud_upload" @click.stop.prevent/>
        </template>
        <template v-slot:append>
          <q-icon name="close" class="cursor-pointer"/>
        </template>
        <template v-slot:hint>
          Upload .fasta files with amino acid sequences
        </template>
      </q-file>
      <PdbViewer v-if="params!.pdbFile" :pdb-file="params!.pdbFile"/>
    </div>
    <div class="col-sm-12 col-md-2"></div>
  </div>
</template>
