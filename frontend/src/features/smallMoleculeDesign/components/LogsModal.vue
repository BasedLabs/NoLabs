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
      required: true
    },
    experimentId: {
      type: String,
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
      if (newShow) {
        this.$q.loading.show({
          spinner: QSpinnerOrbit,
          message: 'Loading logs'
        });
        try {
          this.content = await this.$options.store.logs(this.experimentId);
        }
        finally{
          this.$q.loading.hide();
        }
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
        <div class="row q-col-gutter-xs" v-if="content">
          <div class="col-xs-12 col-md-4">
            Job output
            <q-input
              v-model="content!.output"
              readonly
              filled
              autogrow
            />
          </div>
          <div class="col-xs-12 col-md-4">
            Docking output
            <q-input
              v-model="content!.dockingOutput"
              readonly
              filled
              autogrow
            />
          </div>
          <div class="col-xs-12 col-md-4">
            Errors
            <q-input
              v-model="content!.errors"
              readonly
              filled
              autogrow
            />
          </div>
        </div>
      </q-card-section>
    </q-card>
  </q-dialog>

</template>


