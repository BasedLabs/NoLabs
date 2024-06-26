<script lang="ts">
import {defineComponent, PropType} from 'vue'
import {QSpinnerOrbit} from "quasar";

export default defineComponent({
  name: "JobHeader",
  props: {
    jobName: {
      type: String,
      required: true
    },
    onJobNameChangeSubmit: {
      type: Function as PropType<(newJobName: string) => Promise<void>>,
      required: true
    }
  },
  data() {
    return {
      jobNameData: this.jobName,
      showInferenceForm: false
    }
  },
  methods: {
    changeJobName() {
      this.$q.dialog({
        color: 'info',
        title: 'Prompt',
        message: 'Enter new job name',
        prompt: {
          model: this.jobNameData,
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
        await this.onJobNameChangeSubmit(data);
        this.jobNameData = data;
        this.$q.loading.hide();
      });
    }
  }
})
</script>

<template>
  <q-header :class="$q.dark.isActive ? 'bg-secondary' : 'bg-black'">
    <q-toolbar>
      <q-toolbar-title>{{ jobNameData }}
        <q-btn round
               @click="changeJobName" color="info" size="sm" flat icon="edit"/>
      </q-toolbar-title>
      <slot/>
    </q-toolbar>
  </q-header>
</template>

<style scoped>

</style>