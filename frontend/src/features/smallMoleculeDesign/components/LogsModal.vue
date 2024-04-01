<script lang="ts">
import {defineComponent} from 'vue'
import useSmallMoleculesDesignStore from "../storage";
import {QSpinnerOrbit} from "quasar";

interface Data {
  content: Logs | null;
}

export default defineComponent({
  name: "LogsModal",
  props: {
    visible: {
      type: Boolean,
      default: false,
      required: true
    },
    jobId: {
      type: String,
      default: null,
      required: true
    }
  },
  emits: ['update:visible'],
  store: useSmallMoleculesDesignStore(),
  data(): Data {
    return {
      content: null as Logs | null,
    }
  },
  computed: {
    show: {
      get() {
        return this.visible;
      },
      set(val: boolean) {
        this.$emit('update:visible', val)
      }
    }
  },
  methods: {
    onClose() {
      this.$emit('update:visible', false);
    }
  },
  watch: {
    async visible(newShow, oldShow) {
      debugger;
      if (newShow) {
        this.$q.loading.show({
          spinner: QSpinnerOrbit,
          message: 'Loading logs'
        });
        this.content = await this.$options.store.logs(this.jobId);
        this.$q.loading.hide();
      }
    }
  }
})
</script>

<template>
  <q-dialog
    v-model="show"
  >
    <q-card style="width: 90vh; max-width: 95vw; height: 90vh; max-height: 90vh">
      <q-card-section>
        <div class="text-h6">Logs</div>
      </q-card-section>

      <q-card-section class="q-pt-none">
        <div class="row" v-if="content">
          <div class="col-xs-12 col-md-4">
            <q-input
              v-model="content!.output"
              readonly
              filled
              autogrow
            />
          </div>
          <div class="col-xs-12 col-md-4">
            <q-input
              v-model="content!.dockingOutput"
              readonly
              filled
              autogrow
            />
          </div>
          <div class="col-xs-12 col-md-4">
            <q-input
              v-model="content!.errors"
              readonly
              filled
              autogrow
            />
          </div>
        </div>
      </q-card-section>

      <q-card-actions align="right" class="text-teal q-card-actions-bottom">
        <q-btn flat label="Close" v-close-popup @click="onClose"/>
      </q-card-actions>
    </q-card>
  </q-dialog>

</template>


