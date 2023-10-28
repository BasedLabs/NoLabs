describe('Main page', () => {
    const inputPdbId = '#proteinFileInput';
    const inputAminoAcidSequenceId = '#inputSequence';
    const submitAminoAcidSequenceId = '#submitInference';
    const localisationImageClass = '.localisation-image';
    const foldingSelector = '#viewport div';

    const gets = []
    const logGet = (request) => {
        const start = Date.now()
        request.continue((response) => {
            const duration = Date.now() - start
            gets.push({url: request.url, duration})
        })
    }

    it('Amino Acid lab opened', () => {
        cy.intercept('*', logGet);
        cy.on('uncaught:exception', () => false);
        cy.visit('amino-acid-lab');
    });

    it('Add experiments', () => {
        cy.intercept('*', logGet);
        cy.get('.btn.btn-outline-success.add-experiments-button', {timeout: 90000}).click();
    });

    it('Select experiment', () => {
        cy.intercept('*', logGet);
        cy.get('.bi.bi-check2.btn.btn.btn-outline-success', {timeout: 90000}).click();
        cy.get('#submitInference', {timeout: 90000}).click();
        cy.wait('@all', {timeout: 90000}).its('response.statusCode').should('be.oneOf', [200]);
        cy.get('#mithochondria-list-item', {timeout: 90000}).should('exist');
    });

    it('Submit inference', () => {
        cy.intercept('*', logGet);
        cy.get('#submitInference', {timeout: 90000}).click();
        cy.wait('@all', {timeout: 90000}).its('response.statusCode').should('be.oneOf', [200]);
        cy.get('#mithochondria-list-item', {timeout: 90000}).should('exist');
    });
})