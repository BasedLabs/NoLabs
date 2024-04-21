<script lang="ts">
import {defineComponent, PropType} from 'vue'
import {QSpinnerOrbit, QVueGlobals} from "quasar";

interface ExperimentMetaData {
  id: string;
  name: string;
  date: string;
}

export default defineComponent({
  name: "ExperimentHeader",
  props: {
    metaData: {
      type: Object as PropType<ExperimentMetaData>,
      required: true
    },
    onExperimentNameChangeSubmit: {
      type: Function as PropType<(newExperimentName: string) => Promise<void>>,
      required: true
    }
  },
  data() {
    return {
      experimentMetaData: this.metaData,
      quasar: null as unknown as QVueGlobals,
      showInferenceForm: false
    }
  },
  methods: {
    changeExperimentName() {
      this.$q.dialog({
        color: 'info',
        title: 'Prompt',
        message: 'Enter new localisation name',
        prompt: {
          model: this.experimentMetaData.name,
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
          message: 'Changing localisation name'
        });
        await this.onExperimentNameChangeSubmit(data);
        this.experimentMetaData.name = data;
        this.$q.loading.hide();
      });
    }
  }
})
</script>

<template>
    <q-toolbar>
      <q-item-label class="text-h5 q-pl-md q-pt-md">Experiment name: {{ experimentMetaData.name }}
        <q-btn round
               @click="changeExperimentName" color="info" size="sm" flat icon="edit"/>
      </q-item-label>
    </q-toolbar>
  <q-toolbar>
  <q-item-label class="text-h7 q-pl-md">Experiment id: {{ experimentMetaData.id }}</q-item-label>
  </q-toolbar>

</template>

<style scoped>

</style>
