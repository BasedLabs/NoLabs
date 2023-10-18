<script>
import store from '../storage';
import { api } from '../storage';
import ExperimentsView from './ExperimentsView.vue';
import DrugTargetLabFormView from './DrugTargetLabFormView.vue';
import DrugTargetLabInferenceView from './DrugTargetLabInferenceView.vue';

export default {
    data() {
        return {
            experimentsState: store.state.drugTarget,
            experimentsApi: api.drugTarget
        }
    },
    components: {
        ExperimentsView
    }
}
</script>

<template>
    <ExperimentsView experimentsState="experimentsState" experimentsApi="experimentsApi">
        <template v-slot:labTitle>
            <h4>Drug target lab</h4>
        </template>
        <template v-slot:labForm>
            <DrugTargetLabFormView/>
        </template>
        <template v-slot:lab>
            <DrugTargetLabInferenceView/>
        </template>
    </ExperimentsView >
</template>