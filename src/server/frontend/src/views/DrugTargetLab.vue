<script>
import store from '../storage';
import { api } from '../storage';
import ExperimentsView from './ExperimentsView.vue';
import DrugTargetLabFormView from './DrugTargetLabFormView.vue';
import DrugTargetLabInferenceView from './DrugTargetLabInferenceView.vue';
import DrugTargetNavigator from './DrugTargetNavigator.vue';

export default {
    data() {
        return {
            experimentsState: store.state.drugTarget,
            experimentsApi: api.drugTarget
        }
    },
    components: {
        ExperimentsView,
        DrugTargetLabFormView,
        DrugTargetLabInferenceView,
        DrugTargetNavigator
    }
}
</script>

<template>
    <ExperimentsView :state="experimentsState" :api="experimentsApi">
        <template v-slot:labTitle>
            <h4>Drug target lab</h4>
        </template>
        <template v-slot:experimentMetaData="experimentMetaData" :api="experimentsApi">
            <DrugTargetNavigator :experiment="experimentMetaData.experiment" :api="experimentsApi"/>
        </template>
    </ExperimentsView >
</template>