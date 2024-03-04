<script lang="ts">
import {defineComponent, PropType} from 'vue'
import {QSpinnerOrbit} from "quasar";

export default defineComponent({
  name: "ExperimentHeader",
  props: {
    experimentName: {
      type: String,
      required: true
    },
    onExperimentNameChangeSubmit: {
      type: Function as PropType<(newExperimentName: string) => Promise<void>>,
      required: true
    }
  },
  data() {
    return {
      experimentNameData: this.experimentName,
      showInferenceForm: false
    }
  },
  methods: {
    changeExperimentName() {
      this.$q.dialog({
        color: 'info',
        title: 'Prompt',
        message: 'Enter new experiment name',
        prompt: {
          model: this.experimentNameData,
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
          message: 'Changing experiment name'
        });
        await this.onExperimentNameChangeSubmit(data);
        this.experimentNameData = data;
        this.$q.loading.hide();
      });
    }
  }
})
</script>

<template>
  <q-header :class="$q.dark.isActive ? 'bg-secondary' : 'bg-black'">
    <q-toolbar>
      <q-toolbar-title>{{ experimentNameData }}
        <q-btn round
               @click="changeExperimentName" color="info" size="sm" flat icon="edit"/>
      </q-toolbar-title>
      <slot/>
    </q-toolbar>
  </q-header>
</template>

<style scoped>

</style>